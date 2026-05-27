#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u32 = std::uint32_t;

namespace {

constexpr u64 kMaxN = 216'289'611'853'439'384ULL;
constexpr u64 kMaxOddChunk = 1ull << 23;
constexpr u64 kMinOddChunk = 128;
constexpr std::size_t kStackOddChunk = 4096;

#ifndef NPRIME_CHUNK_FACTOR
#define NPRIME_CHUNK_FACTOR 0.56L
#endif

#ifndef NPRIME_CHUNK_PADDING
#define NPRIME_CHUNK_PADDING 48ULL
#endif

u64 parse_u64(const std::string& value) {
    if (value.empty()) {
        throw std::invalid_argument("empty integer");
    }

    u64 n = 0;
    for (char ch : value) {
        if (ch < '0' || ch > '9') {
            throw std::invalid_argument("invalid decimal integer: " + value);
        }
        const u64 digit = static_cast<u64>(ch - '0');
        if (n > (kMaxN - digit) / 10) {
            throw std::out_of_range("n exceeds supported maximum 216289611853439384");
        }
        n = n * 10 + digit;
    }

    return n;
}

u64 isqrt(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while (r + 1 > r && r + 1 <= x / (r + 1)) {
        ++r;
    }
    while (r > x / r) {
        --r;
    }
    return r;
}

u64 icbrt(u64 x) {
    u64 r = static_cast<u64>(std::cbrt(static_cast<long double>(x)));
    while ((r + 1) <= x / (r + 1) / (r + 1)) {
        ++r;
    }
    while (r > x / r / r) {
        --r;
    }
    return r;
}

u64 iroot4(u64 x) {
    return isqrt(isqrt(x));
}

std::vector<u32> simple_primes(u64 limit) {
    if (limit < 2) {
        return {};
    }
    if (limit > std::numeric_limits<u32>::max()) {
        throw std::overflow_error("sieve limit exceeds uint32 range");
    }

    const u32 n = static_cast<u32>(limit);
    std::vector<std::uint8_t> composite(static_cast<std::size_t>(n) + 1);
    std::vector<u32> primes;
    primes.reserve(static_cast<std::size_t>(limit / std::max<long double>(1.0L, std::log(limit))));

    for (u32 i = 2; i <= n; ++i) {
        if (!composite[i]) {
            primes.push_back(i);
            if (static_cast<u64>(i) * i <= n) {
                for (u64 j = static_cast<u64>(i) * i; j <= n; j += i) {
                    composite[static_cast<std::size_t>(j)] = 1;
                }
            }
        }
    }
    return primes;
}

struct SieveTable {
    u32 limit = 0;
    std::vector<u32> primes;
    std::vector<u32> pi;

    void init(u32 requested) {
        if (limit >= requested) {
            return;
        }

        std::vector<std::uint8_t> composite(static_cast<std::size_t>(requested) + 1);
        primes.clear();
        pi.assign(static_cast<std::size_t>(requested) + 1, 0);

        for (u32 i = 2; i <= requested; ++i) {
            if (!composite[i]) {
                primes.push_back(i);
                if (static_cast<u64>(i) * i <= requested) {
                    for (u64 j = static_cast<u64>(i) * i; j <= requested; j += i) {
                        composite[static_cast<std::size_t>(j)] = 1;
                    }
                }
            }
            pi[i] = static_cast<u32>(primes.size());
        }

        limit = requested;
    }
};

class LehmerCounter {
public:
    LehmerCounter() {
        table_.init(kSieveLimit);
        init_phi_cache();
    }

    u64 count(u64 x) {
        if (x <= table_.limit) {
            return table_.pi[static_cast<std::size_t>(x)];
        }

        const auto found = count_cache_.find(x);
        if (found != count_cache_.end()) {
            return found->second;
        }

        const u64 a = count(iroot4(x));
        const u64 b = count(isqrt(x));
        const u64 c = count(icbrt(x));

        u64 sum = phi(x, a) + ((b + a - 2) * (b - a + 1)) / 2;
        for (u64 i = a; i < b; ++i) {
            const u64 w = x / table_.primes[static_cast<std::size_t>(i)];
            sum -= count(w);

            if (i < c) {
                const u64 bi = count(isqrt(w));
                for (u64 j = i; j < bi; ++j) {
                    sum -= count(w / table_.primes[static_cast<std::size_t>(j)]) - j;
                }
            }
        }

        count_cache_.emplace(x, sum);
        return sum;
    }

    const std::vector<u32>& primes() const {
        return table_.primes;
    }

    u32 limit() const {
        return table_.limit;
    }

private:
    static constexpr u32 kSieveLimit = 1'000'000;
    static constexpr u32 kPhiX = 100'000;
    static constexpr u32 kPhiA = 100;

    SieveTable table_;
    std::vector<u32> phi_cache_;
    std::unordered_map<u64, u64> count_cache_;

    void init_phi_cache() {
        phi_cache_.assign(static_cast<std::size_t>(kPhiA) * kPhiX, 0);
        for (u32 x = 0; x < kPhiX; ++x) {
            phi_cache_[x] = x;
        }
        for (u32 a = 1; a < kPhiA; ++a) {
            const u32 p = table_.primes[static_cast<std::size_t>(a - 1)];
            const std::size_t row = static_cast<std::size_t>(a) * kPhiX;
            const std::size_t prev = static_cast<std::size_t>(a - 1) * kPhiX;
            for (u32 x = 0; x < kPhiX; ++x) {
                phi_cache_[row + x] = phi_cache_[prev + x] - phi_cache_[prev + x / p];
            }
        }
    }

    u64 phi(u64 x, u64 a) {
        if (a == 0) {
            return x;
        }
        if (a < kPhiA && x < kPhiX) {
            return phi_cache_[static_cast<std::size_t>(a) * kPhiX + static_cast<std::size_t>(x)];
        }
        if (x <= table_.limit && a < table_.primes.size()) {
            const u64 next = table_.primes[static_cast<std::size_t>(a)];
            if (x < next * next) {
                return table_.pi[static_cast<std::size_t>(x)] - a + 1;
            }
        }
        return phi(x, a - 1) - phi(x / table_.primes[static_cast<std::size_t>(a - 1)], a - 1);
    }
};

LehmerCounter& lehmer_counter() {
    static LehmerCounter counter;
    return counter;
}

u64 nth_prime_estimate(u64 n) {
    if (n < 6) {
        static constexpr u64 small[] = {0, 2, 3, 5, 7, 11};
        return small[n];
    }

    const long double nn = static_cast<long double>(n);
    const long double l = std::log(nn);
    const long double ll = std::log(l);
    const long double l2 = l * l;
    const long double l3 = l2 * l;
    const long double l4 = l3 * l;
    const long double ll2 = ll * ll;
    const long double ll3 = ll2 * ll;
    const long double ll4 = ll2 * ll2;

    const long double series =
        l + ll - 1.0L +
        (ll - 2.0L) / l -
        (ll2 - 6.0L * ll + 11.0L) / (2.0L * l2) +
        (ll3 - 9.0L * ll2 + 23.0L * ll - 11.0L) / (6.0L * l3) +
        (2.0L * ll4 - 60.0L * ll3 + 570.0L * ll2 - 1920.0L * ll + 2309.0L) / (24.0L * l4);

    return static_cast<u64>(nn * series);
}

u64 initial_radius(u64 n) {
    if (n < 1'000) {
        return 256;
    }

    const long double nn = static_cast<long double>(n);
    const long double l = std::log(nn);
    const long double radius = std::max<long double>(2'000'000.0L, nn / (l * l * l));
    return static_cast<u64>(std::min<long double>(radius, 100'000'000.0L));
}

u64 prime_pi(u64 x) {
    return lehmer_counter().count(x);
}

struct SegmentEntry {
    u32 first_index;
    u32 step;
};

struct SegmentPlan {
    u64 low = 0;
    u64 high = 0;
    u64 odds = 0;
    std::vector<SegmentEntry> entries;
    bool valid = false;
};

const std::vector<SegmentEntry>& segment_plan(const std::vector<u32>& primes, u64 low, u64 high, u64 odds) {
    thread_local SegmentPlan plan;
    if (plan.valid && plan.low == low && plan.high == high && plan.odds == odds) {
        return plan.entries;
    }

    plan.low = low;
    plan.high = high;
    plan.odds = odds;
    plan.entries.clear();

    const u64 root = isqrt(high);
    std::size_t prime_begin = (!primes.empty() && primes.front() == 2) ? 1 : 0;
    const std::size_t prime_end = static_cast<std::size_t>(
        std::upper_bound(primes.begin(), primes.end(), static_cast<u32>(root)) - primes.begin());
    if (prime_end <= prime_begin) {
        plan.valid = true;
        return plan.entries;
    }
    plan.entries.reserve(prime_end - prime_begin);

    for (std::size_t prime_idx = prime_begin; prime_idx < prime_end; ++prime_idx) {
        const u32 p = primes[prime_idx];
        u64 first = static_cast<u64>(p) * p;
        if (first < low) {
            const u64 x = low + p - 1;
            first = (x / p) * static_cast<u64>(p);
        }
        if ((first & 1u) == 0) {
            first += p;
        }
        if (first <= high) {
            plan.entries.push_back({
                static_cast<u32>((first - low) / 2),
                p
            });
        }
    }

    plan.valid = true;
    return plan.entries;
}

u64 nth_prime_after(const std::vector<u32>& primes, u64 start, u64 offset) {
    u64 seen = 0;
    u64 chunk_low = std::max<u64>(start + 1, 3);
    if (offset == 0) {
        throw std::logic_error("forward prime offset must be positive");
    }
    if (start < 2) {
        ++seen;
        if (seen == offset) {
            return 2;
        }
    }
    if ((chunk_low & 1u) == 0) {
        ++chunk_low;
    }

    for (;;) {
        const long double gap = std::max<long double>(8.0L, std::log(static_cast<long double>(chunk_low)));
        const u64 remaining = offset - seen;
        const u64 wanted =
            static_cast<u64>(static_cast<long double>(remaining) * gap * NPRIME_CHUNK_FACTOR) +
            NPRIME_CHUNK_PADDING;
        const u64 chunk_odds = std::clamp<u64>(wanted, kMinOddChunk, kMaxOddChunk);
        const u64 chunk_high = chunk_low + 2 * (chunk_odds - 1);
        std::array<std::uint8_t, kStackOddChunk> stack_composite;
        thread_local std::vector<std::uint8_t> composite;
        std::uint8_t* composite_data = stack_composite.data();
        if (chunk_odds <= kStackOddChunk) {
            std::fill_n(composite_data, static_cast<std::size_t>(chunk_odds), 0);
        } else {
            composite.assign(static_cast<std::size_t>(chunk_odds), 0);
            composite_data = composite.data();
        }

        for (const SegmentEntry& entry : segment_plan(primes, chunk_low, chunk_high, chunk_odds)) {
            for (u32 index = entry.first_index; index < chunk_odds; index += entry.step) {
                composite_data[index] = 1;
            }
        }

        for (u64 i = 0; i < chunk_odds; ++i) {
            if (!composite_data[static_cast<std::size_t>(i)]) {
                ++seen;
                if (seen == offset) {
                    return chunk_low + 2 * i;
                }
            }
        }

        if (chunk_high > std::numeric_limits<u64>::max() - 2) {
            break;
        }
        chunk_low = chunk_high + 2;
    }

    throw std::overflow_error("nth prime exceeds uint64 range");
}

u64 nth_prime_before_or_at(const std::vector<u32>& primes, u64 start, u64 skip) {
    u64 seen = 0;
    u64 chunk_high = start;
    if (chunk_high >= 2 && chunk_high < 3) {
        return 2;
    }
    if ((chunk_high & 1u) == 0) {
        --chunk_high;
    }

    while (chunk_high >= 3) {
        const long double gap = std::max<long double>(8.0L, std::log(static_cast<long double>(chunk_high)));
        const u64 remaining = skip - std::min(seen, skip);
        const u64 wanted =
            static_cast<u64>(static_cast<long double>(remaining) * gap * NPRIME_CHUNK_FACTOR) +
            NPRIME_CHUNK_PADDING;
        const u64 chunk_odds = std::min<u64>(
            ((chunk_high - 3) / 2) + 1,
            std::clamp<u64>(wanted, kMinOddChunk, kMaxOddChunk));
        const u64 chunk_low = chunk_high - 2 * (chunk_odds - 1);
        std::array<std::uint8_t, kStackOddChunk> stack_composite;
        thread_local std::vector<std::uint8_t> composite;
        std::uint8_t* composite_data = stack_composite.data();
        if (chunk_odds <= kStackOddChunk) {
            std::fill_n(composite_data, static_cast<std::size_t>(chunk_odds), 0);
        } else {
            composite.assign(static_cast<std::size_t>(chunk_odds), 0);
            composite_data = composite.data();
        }

        for (const SegmentEntry& entry : segment_plan(primes, chunk_low, chunk_high, chunk_odds)) {
            for (u32 index = entry.first_index; index < chunk_odds; index += entry.step) {
                composite_data[index] = 1;
            }
        }

        for (u64 i = chunk_odds; i-- > 0;) {
            if (!composite_data[static_cast<std::size_t>(i)]) {
                if (seen == skip) {
                    return chunk_low + 2 * i;
                }
                ++seen;
            }
        }

        if (chunk_low <= 3) {
            break;
        }
        chunk_high = chunk_low - 2;
    }

    if (seen == skip) {
        return 2;
    }

    throw std::logic_error("target prime was not found inside bracket");
}

u64 nth_prime_counting(u64 n) {
    struct SearchState {
        u64 n = 0;
        u64 estimate = 0;
        u64 count_at_estimate = 0;
        bool table_primes_ok = false;
        bool valid = false;
    };

    thread_local SearchState cached_state;

    u64 estimate = 0;
    u64 count_at_estimate = 0;
    bool table_primes_ok = false;
    if (cached_state.valid && cached_state.n == n) {
        estimate = cached_state.estimate;
        count_at_estimate = cached_state.count_at_estimate;
        table_primes_ok = cached_state.table_primes_ok;
    } else {
        estimate = nth_prime_estimate(n);
        count_at_estimate = prime_pi(estimate);

        const u64 delta = count_at_estimate > n ? count_at_estimate - n : n - count_at_estimate;
        if (delta > 256) {
            const long double signed_delta = count_at_estimate >= n
                ? static_cast<long double>(count_at_estimate - n)
                : -static_cast<long double>(n - count_at_estimate);
            long double adjusted = static_cast<long double>(estimate) - signed_delta * std::log(static_cast<long double>(estimate));
            if (adjusted < 2.0L) {
                adjusted = 2.0L;
            }

            const u64 adjusted_estimate = static_cast<u64>(adjusted);
            if (adjusted_estimate != estimate) {
                const u64 adjusted_count = prime_pi(adjusted_estimate);
                const u64 adjusted_delta = adjusted_count > n ? adjusted_count - n : n - adjusted_count;
                if (adjusted_delta < delta) {
                    estimate = adjusted_estimate;
                    count_at_estimate = adjusted_count;
                }
            }
        }

        const u64 needed_limit = isqrt(estimate + initial_radius(n)) + 1;
        table_primes_ok = needed_limit <= lehmer_counter().limit();
        cached_state = {n, estimate, count_at_estimate, table_primes_ok, true};
    }

    const std::vector<u32>* primes = &lehmer_counter().primes();
    std::vector<u32> generated_primes;
    if (!table_primes_ok) {
        const u64 needed_limit = isqrt(estimate + initial_radius(n)) + 1;
        generated_primes = simple_primes(needed_limit);
        primes = &generated_primes;
    }

    if (count_at_estimate < n) {
        return nth_prime_after(*primes, estimate, n - count_at_estimate);
    }
    return nth_prime_before_or_at(*primes, estimate, count_at_estimate - n);
}

void trace_nth_prime(u64 n) {
    if (n == 0) {
        std::cout << "n=0 result=0\n";
        return;
    }

    u64 estimate = nth_prime_estimate(n);
    u64 count_at_estimate = prime_pi(estimate);
    const u64 initial_delta = count_at_estimate > n ? count_at_estimate - n : n - count_at_estimate;
    std::cout << "initial_estimate=" << estimate
              << " pi=" << count_at_estimate
              << " delta=" << (count_at_estimate >= n ? "+" : "-") << initial_delta << '\n';

    if (initial_delta > 256) {
        const long double signed_delta = count_at_estimate >= n
            ? static_cast<long double>(count_at_estimate - n)
            : -static_cast<long double>(n - count_at_estimate);
        long double adjusted = static_cast<long double>(estimate) - signed_delta * std::log(static_cast<long double>(estimate));
        if (adjusted < 2.0L) {
            adjusted = 2.0L;
        }

        const u64 adjusted_estimate = static_cast<u64>(adjusted);
        const u64 adjusted_count = prime_pi(adjusted_estimate);
        const u64 adjusted_delta = adjusted_count > n ? adjusted_count - n : n - adjusted_count;
        std::cout << "adjusted_estimate=" << adjusted_estimate
                  << " pi=" << adjusted_count
                  << " delta=" << (adjusted_count >= n ? "+" : "-") << adjusted_delta << '\n';
    }

    const auto start = std::chrono::steady_clock::now();
    const u64 result = nth_prime_counting(n);
    const auto stop = std::chrono::steady_clock::now();
    const auto elapsed =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(stop - start);
    std::cout << "result=" << result << " elapsed_ms=" << elapsed.count() << '\n';
}

u64 nth_prime(u64 n) {
    if (n == 0) {
        return 0;
    }
    return nth_prime_counting(n);
}

bool self_test() {
    const std::vector<std::pair<u64, u64>> cases = {
        {0, 0},
        {1, 2},
        {2, 3},
        {3, 5},
        {4, 7},
        {5, 11},
        {6, 13},
        {10, 29},
        {100, 541},
        {1'000, 7'919},
        {10'000, 104'729},
        {100'000, 1'299'709},
        {1'000'000, 15'485'863},
        {1'000'000'000, 22'801'763'489ULL},
        {10'000'000'000ULL, 252'097'800'623ULL},
    };

    for (const auto& test : cases) {
        const u64 got = nth_prime(test.first);
        if (got != test.second) {
            std::cerr << "self-test failed for n=" << test.first
                      << ": got " << got << ", expected " << test.second << '\n';
            return false;
        }
    }
    return true;
}

void usage(const char* argv0) {
    std::cerr << "usage:\n"
              << "  " << argv0 << " N\n"
              << "  " << argv0 << " --bench N [iterations]\n"
              << "  " << argv0 << " --trace N\n"
              << "  " << argv0 << " --self-test\n";
}

template <typename T>
void benchmark_barrier(T& value) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : "+r"(value) : : "memory");
#else
    volatile T* sink = &value;
    (void)sink;
#endif
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            if (!self_test()) {
                return 1;
            }
            std::cout << "ok\n";
            return 0;
        }

        if (argc >= 3 && std::string(argv[1]) == "--bench") {
            const u64 n = parse_u64(argv[2]);
            const u64 iterations = argc >= 4 ? parse_u64(argv[3]) : 10;
            if (iterations == 0) {
                throw std::invalid_argument("iterations must be positive");
            }

            u64 result = 0;
            const auto start = std::chrono::steady_clock::now();
            for (u64 i = 0; i < iterations; ++i) {
                u64 query = n;
                benchmark_barrier(query);
                result = nth_prime(query);
                benchmark_barrier(result);
            }
            const auto stop = std::chrono::steady_clock::now();
            const auto elapsed =
                std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(stop - start);

            std::cout << result << '\n'
                      << "iterations: " << iterations << '\n'
                      << "total_ms: " << std::fixed << std::setprecision(3) << elapsed.count() << '\n'
                      << "avg_ms: " << std::fixed << std::setprecision(6)
                      << elapsed.count() / static_cast<double>(iterations) << '\n';
            return 0;
        }

        if (argc == 3 && std::string(argv[1]) == "--trace") {
            trace_nth_prime(parse_u64(argv[2]));
            return 0;
        }

        if (argc != 2) {
            usage(argv[0]);
            return 2;
        }

        std::cout << nth_prime(parse_u64(argv[1])) << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
