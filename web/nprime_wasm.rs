#![allow(static_mut_refs)]

use std::collections::HashMap;

type U64 = u64;
type U32 = u32;

const SIEVE_LIMIT: U32 = 1_000_000;
const PHI_X: U32 = 100_000;
const PHI_A: U32 = 100;
const MAX_ODD_CHUNK: U64 = 1 << 23;
const MIN_ODD_CHUNK: U64 = 128;
const STACK_ODD_CHUNK: usize = 4096;
const CHUNK_FACTOR: f64 = 0.56;
const CHUNK_PADDING: U64 = 48;

struct SieveTable {
    limit: U32,
    primes: Vec<U32>,
    pi: Vec<U32>,
}

impl SieveTable {
    fn new() -> Self {
        Self {
            limit: 0,
            primes: Vec::new(),
            pi: Vec::new(),
        }
    }

    fn init(&mut self, requested: U32) {
        if self.limit >= requested {
            return;
        }

        let mut composite = vec![0u8; requested as usize + 1];
        self.primes.clear();
        self.pi.clear();
        self.pi.resize(requested as usize + 1, 0);

        for i in 2..=requested {
            if composite[i as usize] == 0 {
                self.primes.push(i);
                let square = i as U64 * i as U64;
                if square <= requested as U64 {
                    let mut j = square;
                    while j <= requested as U64 {
                        composite[j as usize] = 1;
                        j += i as U64;
                    }
                }
            }
            self.pi[i as usize] = self.primes.len() as U32;
        }

        self.limit = requested;
    }
}

struct LehmerCounter {
    table: SieveTable,
    phi_cache: Vec<U32>,
    count_cache: HashMap<U64, U64>,
}

impl LehmerCounter {
    fn new() -> Self {
        let mut table = SieveTable::new();
        table.init(SIEVE_LIMIT);
        let mut counter = Self {
            table,
            phi_cache: Vec::new(),
            count_cache: HashMap::new(),
        };
        counter.init_phi_cache();
        counter
    }

    fn init_phi_cache(&mut self) {
        self.phi_cache.clear();
        self.phi_cache.resize(PHI_A as usize * PHI_X as usize, 0);

        for x in 0..PHI_X {
            self.phi_cache[x as usize] = x;
        }
        for a in 1..PHI_A {
            let p = self.table.primes[(a - 1) as usize];
            let row = a as usize * PHI_X as usize;
            let prev = (a - 1) as usize * PHI_X as usize;
            for x in 0..PHI_X {
                self.phi_cache[row + x as usize] =
                    self.phi_cache[prev + x as usize] - self.phi_cache[prev + (x / p) as usize];
            }
        }
    }

    fn count(&mut self, x: U64) -> U64 {
        if x <= self.table.limit as U64 {
            return self.table.pi[x as usize] as U64;
        }
        if let Some(value) = self.count_cache.get(&x) {
            return *value;
        }

        let a = self.count(iroot4(x));
        let b = self.count(isqrt(x));
        let c = self.count(icbrt(x));

        let mut sum = self.phi(x, a) + ((b + a - 2) * (b - a + 1)) / 2;
        for i in a..b {
            let p_i = self.table.primes[i as usize] as U64;
            let w = x / p_i;
            sum -= self.count(w);

            if i < c {
                let bi = self.count(isqrt(w));
                for j in i..bi {
                    let p_j = self.table.primes[j as usize] as U64;
                    sum -= self.count(w / p_j) - j;
                }
            }
        }

        self.count_cache.insert(x, sum);
        sum
    }

    fn phi(&self, x: U64, a: U64) -> U64 {
        if a == 0 {
            return x;
        }
        if a < PHI_A as U64 && x < PHI_X as U64 {
            return self.phi_cache[a as usize * PHI_X as usize + x as usize] as U64;
        }
        if x <= self.table.limit as U64 && (a as usize) < self.table.primes.len() {
            let next = self.table.primes[a as usize] as U64;
            if x < next * next {
                return self.table.pi[x as usize] as U64 - a + 1;
            }
        }
        self.phi(x, a - 1) - self.phi(x / self.table.primes[(a - 1) as usize] as U64, a - 1)
    }
}

#[derive(Clone, Copy)]
struct SegmentEntry {
    first_index: U32,
    step: U32,
}

struct SegmentPlan {
    low: U64,
    high: U64,
    odds: U64,
    entries: Vec<SegmentEntry>,
    valid: bool,
}

impl SegmentPlan {
    fn new() -> Self {
        Self {
            low: 0,
            high: 0,
            odds: 0,
            entries: Vec::new(),
            valid: false,
        }
    }
}

#[derive(Clone, Copy)]
struct SearchState {
    n: U64,
    estimate: U64,
    count_at_estimate: U64,
    table_primes_ok: bool,
    valid: bool,
}

impl SearchState {
    const fn empty() -> Self {
        Self {
            n: 0,
            estimate: 0,
            count_at_estimate: 0,
            table_primes_ok: false,
            valid: false,
        }
    }
}

static mut COUNTER: Option<LehmerCounter> = None;
static mut SEGMENT_PLAN: Option<SegmentPlan> = None;
static mut SEARCH_STATE: SearchState = SearchState::empty();

fn counter() -> &'static mut LehmerCounter {
    unsafe {
        if COUNTER.is_none() {
            COUNTER = Some(LehmerCounter::new());
        }
        COUNTER.as_mut().unwrap()
    }
}

fn isqrt(x: U64) -> U64 {
    let mut r = (x as f64).sqrt() as U64;
    while r + 1 > r && r + 1 <= x / (r + 1) {
        r += 1;
    }
    while r > x / r {
        r -= 1;
    }
    r
}

fn icbrt(x: U64) -> U64 {
    let mut r = (x as f64).cbrt() as U64;
    while r + 1 <= x / (r + 1) / (r + 1) {
        r += 1;
    }
    while r > x / r / r {
        r -= 1;
    }
    r
}

fn iroot4(x: U64) -> U64 {
    isqrt(isqrt(x))
}

fn simple_primes(limit: U64) -> Vec<U32> {
    if limit < 2 {
        return Vec::new();
    }

    let n = limit as U32;
    let mut composite = vec![0u8; n as usize + 1];
    let reserve = (limit as f64 / (limit as f64).ln().max(1.0)).max(8.0) as usize;
    let mut primes = Vec::with_capacity(reserve);

    for i in 2..=n {
        if composite[i as usize] == 0 {
            primes.push(i);
            let square = i as U64 * i as U64;
            if square <= n as U64 {
                let mut j = square;
                while j <= n as U64 {
                    composite[j as usize] = 1;
                    j += i as U64;
                }
            }
        }
    }
    primes
}

fn nth_prime_estimate(n: U64) -> U64 {
    match n {
        0 => 0,
        1 => 2,
        2 => 3,
        3 => 5,
        4 => 7,
        5 => 11,
        _ => {
            let nn = n as f64;
            let l = nn.ln();
            let ll = l.ln();
            let l2 = l * l;
            let l3 = l2 * l;
            let l4 = l3 * l;
            let ll2 = ll * ll;
            let ll3 = ll2 * ll;
            let ll4 = ll2 * ll2;

            let series = l
                + ll
                - 1.0
                + (ll - 2.0) / l
                - (ll2 - 6.0 * ll + 11.0) / (2.0 * l2)
                + (ll3 - 9.0 * ll2 + 23.0 * ll - 11.0) / (6.0 * l3)
                + (2.0 * ll4 - 60.0 * ll3 + 570.0 * ll2 - 1920.0 * ll + 2309.0)
                    / (24.0 * l4);

            (nn * series) as U64
        }
    }
}

fn initial_radius(n: U64) -> U64 {
    if n < 1_000 {
        return 256;
    }
    let nn = n as f64;
    let l = nn.ln();
    (2_000_000.0f64.max(nn / (l * l * l))).min(100_000_000.0) as U64
}

fn prime_pi(x: U64) -> U64 {
    counter().count(x)
}

fn segment_plan(primes: &[U32], low: U64, high: U64, odds: U64) -> &'static [SegmentEntry] {
    unsafe {
        if SEGMENT_PLAN.is_none() {
            SEGMENT_PLAN = Some(SegmentPlan::new());
        }
        let plan = SEGMENT_PLAN.as_mut().unwrap();
        if plan.valid && plan.low == low && plan.high == high && plan.odds == odds {
            return &plan.entries;
        }

        plan.low = low;
        plan.high = high;
        plan.odds = odds;
        plan.entries.clear();

        let root = isqrt(high);
        let prime_begin = if !primes.is_empty() && primes[0] == 2 { 1 } else { 0 };
        let mut prime_end = prime_begin;
        while prime_end < primes.len() && primes[prime_end] as U64 <= root {
            prime_end += 1;
        }

        if prime_end > prime_begin {
            plan.entries.reserve(prime_end - prime_begin);
            for &p32 in &primes[prime_begin..prime_end] {
                let p = p32 as U64;
                let mut first = p * p;
                if first < low {
                    let x = low + p - 1;
                    first = (x / p) * p;
                }
                if first & 1 == 0 {
                    first += p;
                }
                if first <= high {
                    plan.entries.push(SegmentEntry {
                        first_index: ((first - low) / 2) as U32,
                        step: p32,
                    });
                }
            }
        }

        plan.valid = true;
        &plan.entries
    }
}

fn nth_prime_after(primes: &[U32], start: U64, offset: U64) -> U64 {
    let mut seen = 0;
    let mut chunk_low = (start + 1).max(3);
    if start < 2 {
        seen += 1;
        if seen == offset {
            return 2;
        }
    }
    if chunk_low & 1 == 0 {
        chunk_low += 1;
    }

    loop {
        let gap = (chunk_low as f64).ln().max(8.0);
        let remaining = offset - seen;
        let wanted = (remaining as f64 * gap * CHUNK_FACTOR) as U64 + CHUNK_PADDING;
        let chunk_odds = wanted.clamp(MIN_ODD_CHUNK, MAX_ODD_CHUNK);
        let chunk_high = chunk_low + 2 * (chunk_odds - 1);

        let mut stack = [0u8; STACK_ODD_CHUNK];
        let mut heap = Vec::new();
        let composite: &mut [u8] = if chunk_odds as usize <= STACK_ODD_CHUNK {
            &mut stack[..chunk_odds as usize]
        } else {
            heap.resize(chunk_odds as usize, 0);
            heap.as_mut_slice()
        };

        for entry in segment_plan(primes, chunk_low, chunk_high, chunk_odds) {
            let mut index = entry.first_index;
            while (index as U64) < chunk_odds {
                composite[index as usize] = 1;
                index += entry.step;
            }
        }

        for i in 0..chunk_odds {
            if composite[i as usize] == 0 {
                seen += 1;
                if seen == offset {
                    return chunk_low + 2 * i;
                }
            }
        }

        chunk_low = chunk_high + 2;
    }
}

fn nth_prime_before_or_at(primes: &[U32], start: U64, skip: U64) -> U64 {
    let mut seen = 0;
    let mut chunk_high = start;
    if chunk_high >= 2 && chunk_high < 3 {
        return 2;
    }
    if chunk_high & 1 == 0 {
        chunk_high -= 1;
    }

    while chunk_high >= 3 {
        let gap = (chunk_high as f64).ln().max(8.0);
        let remaining = skip - seen.min(skip);
        let wanted = (remaining as f64 * gap * CHUNK_FACTOR) as U64 + CHUNK_PADDING;
        let cap = ((chunk_high - 3) / 2) + 1;
        let chunk_odds = cap.min(wanted.clamp(MIN_ODD_CHUNK, MAX_ODD_CHUNK));
        let chunk_low = chunk_high - 2 * (chunk_odds - 1);

        let mut stack = [0u8; STACK_ODD_CHUNK];
        let mut heap = Vec::new();
        let composite: &mut [u8] = if chunk_odds as usize <= STACK_ODD_CHUNK {
            &mut stack[..chunk_odds as usize]
        } else {
            heap.resize(chunk_odds as usize, 0);
            heap.as_mut_slice()
        };

        for entry in segment_plan(primes, chunk_low, chunk_high, chunk_odds) {
            let mut index = entry.first_index;
            while (index as U64) < chunk_odds {
                composite[index as usize] = 1;
                index += entry.step;
            }
        }

        let mut i = chunk_odds;
        while i > 0 {
            i -= 1;
            if composite[i as usize] == 0 {
                if seen == skip {
                    return chunk_low + 2 * i;
                }
                seen += 1;
            }
        }

        if chunk_low <= 3 {
            break;
        }
        chunk_high = chunk_low - 2;
    }

    2
}

fn nth_prime_counting(n: U64) -> U64 {
    let (estimate, count_at_estimate, table_primes_ok) = unsafe {
        if SEARCH_STATE.valid && SEARCH_STATE.n == n {
            (
                SEARCH_STATE.estimate,
                SEARCH_STATE.count_at_estimate,
                SEARCH_STATE.table_primes_ok,
            )
        } else {
            let mut estimate = nth_prime_estimate(n);
            let mut count_at_estimate = prime_pi(estimate);
            let delta = count_at_estimate.abs_diff(n);

            if delta > 256 {
                let signed_delta = if count_at_estimate >= n {
                    (count_at_estimate - n) as f64
                } else {
                    -((n - count_at_estimate) as f64)
                };
                let mut adjusted = estimate as f64 - signed_delta * (estimate as f64).ln();
                if adjusted < 2.0 {
                    adjusted = 2.0;
                }

                let adjusted_estimate = adjusted as U64;
                if adjusted_estimate != estimate {
                    let adjusted_count = prime_pi(adjusted_estimate);
                    if adjusted_count.abs_diff(n) < delta {
                        estimate = adjusted_estimate;
                        count_at_estimate = adjusted_count;
                    }
                }
            }

            let needed_limit = isqrt(estimate + initial_radius(n)) + 1;
            let table_primes_ok = needed_limit <= counter().table.limit as U64;
            SEARCH_STATE = SearchState {
                n,
                estimate,
                count_at_estimate,
                table_primes_ok,
                valid: true,
            };
            (estimate, count_at_estimate, table_primes_ok)
        }
    };

    let generated_primes;
    let primes: &[U32] = if table_primes_ok {
        &counter().table.primes
    } else {
        let needed_limit = isqrt(estimate + initial_radius(n)) + 1;
        generated_primes = simple_primes(needed_limit);
        &generated_primes
    };

    if count_at_estimate < n {
        nth_prime_after(primes, estimate, n - count_at_estimate)
    } else {
        nth_prime_before_or_at(primes, estimate, count_at_estimate - n)
    }
}

fn nth_prime_with_runtime_state(n: U64) -> U64 {
    if n == 0 {
        0
    } else {
        nth_prime_counting(n)
    }
}

#[no_mangle]
pub extern "C" fn reset_state() {
    unsafe {
        COUNTER = None;
        SEGMENT_PLAN = None;
        SEARCH_STATE = SearchState::empty();
    }
}

#[no_mangle]
pub extern "C" fn nth_prime(n: U64) -> U64 {
    reset_state();
    let result = nth_prime_with_runtime_state(n);
    reset_state();
    result
}

#[no_mangle]
pub extern "C" fn nth_prime_isolated(n: U64) -> U64 {
    nth_prime(n)
}

#[no_mangle]
pub extern "C" fn estimate_nth_prime(n: U64) -> U64 {
    nth_prime_estimate(n)
}

#[no_mangle]
pub extern "C" fn prime_pi_exact(x: U64) -> U64 {
    reset_state();
    let result = prime_pi(x);
    reset_state();
    result
}

#[no_mangle]
pub extern "C" fn self_test_small() -> U32 {
    let cases: &[(U64, U64)] = &[
        (0, 0),
        (1, 2),
        (2, 3),
        (3, 5),
        (4, 7),
        (5, 11),
        (6, 13),
        (10, 29),
        (100, 541),
        (1_000, 7_919),
        (10_000, 104_729),
        (100_000, 1_299_709),
        (1_000_000, 15_485_863),
    ];

    for &(n, expected) in cases {
        if nth_prime(n) != expected {
            return 0;
        }
    }
    1
}

#[no_mangle]
pub extern "C" fn algorithm_revision() -> U32 {
    1
}
