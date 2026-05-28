#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u32 = std::uint32_t;
using u16 = std::uint16_t;

namespace {

constexpr u64 kMaxN = 216'289'611'853'439'384ULL;
constexpr u64 kMaxOddChunk = 1ull << 23;
constexpr u64 kMinOddChunk = 128;
constexpr std::size_t kStackOddChunk = 65536;
constexpr u32 kWheelPrimeCount = 6;
constexpr u32 kRankBlockOdds = 64;

#ifndef NPRIME_CHUNK_FACTOR
#define NPRIME_CHUNK_FACTOR 0.56L
#endif

#ifndef NPRIME_CHUNK_PADDING
#define NPRIME_CHUNK_PADDING 48ULL
#endif

#ifndef NPRIME_PHI_X
#define NPRIME_PHI_X 100000U
#endif

#ifndef NPRIME_PHI_A
#define NPRIME_PHI_A 100U
#endif

#ifndef NPRIME_SIEVE_LIMIT
#define NPRIME_SIEVE_LIMIT 1000000U
#endif

#ifndef NPRIME_DIRECT_WHEEL_A
#define NPRIME_DIRECT_WHEEL_A 10
#endif

#ifndef NPRIME_TOP_CHUNK_BIG
#define NPRIME_TOP_CHUNK_BIG 32ULL
#endif

#ifndef NPRIME_TOP_CHUNK_MID
#define NPRIME_TOP_CHUNK_MID 32ULL
#endif

#ifndef NPRIME_TOP_CHUNK_SMALL
#define NPRIME_TOP_CHUNK_SMALL 64ULL
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

u32 parse_u32_arg(const std::string& value, const char* name) {
    const u64 parsed = parse_u64(value);
    if (parsed > std::numeric_limits<u32>::max()) {
        throw std::out_of_range(std::string(name) + " exceeds uint32 range");
    }
    return static_cast<u32>(parsed);
}

u64 isqrt(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<double>(x)));
    while (r + 1 > r && r + 1 <= x / (r + 1)) {
        ++r;
    }
    while (r > x / r) {
        --r;
    }
    return r;
}

u64 icbrt(u64 x) {
    u64 r = static_cast<u64>(std::cbrt(static_cast<double>(x)));
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

u32 popcount64(u64 value) {
#if defined(__GNUC__) || defined(__clang__)
    return static_cast<u32>(__builtin_popcountll(value));
#else
    u32 count = 0;
    while (value != 0) {
        value &= value - 1;
        ++count;
    }
    return count;
#endif
}

double elapsed_ms(
    std::chrono::steady_clock::time_point start,
    std::chrono::steady_clock::time_point stop) {
    return std::chrono::duration<double, std::milli>(stop - start).count();
}

u64 elapsed_ns(
    std::chrono::steady_clock::time_point start,
    std::chrono::steady_clock::time_point stop) {
    return static_cast<u64>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count());
}

std::vector<u32> simple_primes(u64 limit) {
    if (limit < 2) {
        return {};
    }
    if (limit > std::numeric_limits<u32>::max()) {
        throw std::overflow_error("sieve limit exceeds uint32 range");
    }

    const u32 n = static_cast<u32>(limit);
    std::vector<std::uint8_t> composite(static_cast<std::size_t>(n / 2) + 1);
    std::vector<u32> primes;
    primes.reserve(static_cast<std::size_t>(limit / std::max<long double>(1.0L, std::log(limit))));
    primes.push_back(2);

    const u32 root = static_cast<u32>(isqrt(n));
    for (u32 i = 3; i <= root; i += 2) {
        if (!composite[i / 2]) {
            for (u64 j = static_cast<u64>(i) * i; j <= n; j += 2ull * i) {
                composite[static_cast<std::size_t>(j / 2)] = 1;
            }
        }
    }

    for (u32 i = 3; i <= n; i += 2) {
        if (!composite[i / 2]) {
            primes.push_back(i);
        }
    }
    return primes;
}

struct SieveBuildTimings {
    double alloc_ms = 0.0;
    double mark_ms = 0.0;
    double pi_ms = 0.0;
    double rank_ms = 0.0;
};

struct SieveTable {
    u32 limit = 0;
    u32 rank_limit = 0;
    std::vector<u32> primes;
    std::vector<u16> pi;
    std::vector<u64> rank_composite;
    std::vector<u32> rank_pi_block;

    void init(u32 requested, u32 rank_requested, SieveBuildTimings* timings = nullptr) {
        rank_requested = std::max(rank_requested, requested);
        if (limit >= requested && rank_limit >= rank_requested) {
            return;
        }

        const u32 sieve_requested = std::max(requested, rank_requested);
        std::chrono::steady_clock::time_point stage_start;
        if (timings != nullptr) {
            stage_start = std::chrono::steady_clock::now();
        }
        std::vector<std::uint8_t> composite(static_cast<std::size_t>(sieve_requested / 2) + 1);
        primes.clear();
        primes.reserve(
            static_cast<std::size_t>(
                1.3L * static_cast<long double>(requested) /
                std::max<long double>(1.0L, std::log(static_cast<long double>(std::max<u32>(requested, 3))))) +
            32);
        pi.assign(static_cast<std::size_t>(requested) + 1, 0);
        rank_composite.clear();
        rank_pi_block.clear();
        if (timings != nullptr) {
            const auto stage_stop = std::chrono::steady_clock::now();
            timings->alloc_ms += elapsed_ms(stage_start, stage_stop);
            stage_start = stage_stop;
        }
        if (sieve_requested < 2) {
            limit = requested;
            rank_limit = rank_requested;
            return;
        }

        const u32 root = static_cast<u32>(isqrt(sieve_requested));
        for (u32 i = 3; i <= root; i += 2) {
            if (!composite[i / 2]) {
                for (u64 j = static_cast<u64>(i) * i; j <= sieve_requested; j += 2ull * i) {
                    composite[static_cast<std::size_t>(j / 2)] = 1;
                }
            }
        }
        if (timings != nullptr) {
            const auto stage_stop = std::chrono::steady_clock::now();
            timings->mark_ms += elapsed_ms(stage_start, stage_stop);
            stage_start = stage_stop;
        }

        primes.push_back(2);
        u16 count = 1;
        if (requested >= 2) {
            pi[2] = count;
        }
        for (u32 i = 3; i <= requested; i += 2) {
            pi[i - 1] = count;
            if (!composite[i / 2]) {
                primes.push_back(i);
                if (count == std::numeric_limits<u16>::max()) {
                    throw std::overflow_error("pi table count exceeds uint16 range");
                }
                ++count;
            }
            pi[i] = count;
        }
        if ((requested & 1u) == 0) {
            pi[requested] = count;
        }
        if (timings != nullptr) {
            const auto stage_stop = std::chrono::steady_clock::now();
            timings->pi_ms += elapsed_ms(stage_start, stage_stop);
            stage_start = stage_stop;
        }

        limit = requested;
        rank_limit = rank_requested;
        if (rank_requested > requested) {
            const u32 odd_count = rank_requested >= 3 ? ((rank_requested - 3) / 2) + 1 : 0;
            const u32 blocks = (odd_count + kRankBlockOdds - 1) / kRankBlockOdds;
            rank_composite.assign(blocks, 0);
            rank_pi_block.assign(static_cast<std::size_t>(blocks) + 1, 0);

            u32 count = 0;
            for (u32 block = 0; block < blocks; ++block) {
                rank_pi_block[block] = count;
                const u32 begin = block * kRankBlockOdds;
                const u32 valid = std::min<u32>(kRankBlockOdds, odd_count - begin);
                u64 bits = 0;
                for (u32 bit = 0; bit < valid; ++bit) {
                    const u32 n = 3 + 2 * (begin + bit);
                    if (composite[n / 2]) {
                        bits |= 1ull << bit;
                    }
                }
                rank_composite[block] = bits;
                count += valid - popcount64(bits);
            }
            rank_pi_block[blocks] = count;
        }
        if (timings != nullptr) {
            timings->rank_ms += elapsed_ms(stage_start, std::chrono::steady_clock::now());
        }
    }

    u32 count(u64 x) const {
        return pi[static_cast<std::size_t>(x)];
    }

    u32 rank_count(u64 x) const {
        if (x <= limit) {
            return count(x);
        }
        const u32 value = static_cast<u32>(x);
        const u32 odd_count = value >= 3 ? ((value - 3) / 2) + 1 : 0;
        const u32 block = odd_count / kRankBlockOdds;
        const u32 rem = odd_count % kRankBlockOdds;
        u32 total = 1 + rank_pi_block[block];
        if (rem != 0) {
            const u64 mask = (1ull << rem) - 1;
            total += rem - popcount64(rank_composite[block] & mask);
        }
        return total;
    }
};

struct CountConfig {
    u32 sieve_limit;
    u32 rank_limit;
    u32 phi_x;
    u32 phi_a;
    u32 threads;
    u32 a_boost;
    u32 phi_split_slack;

    bool operator==(const CountConfig& other) const {
        return sieve_limit == other.sieve_limit && rank_limit == other.rank_limit && phi_x == other.phi_x &&
            phi_a == other.phi_a && threads == other.threads &&
            a_boost == other.a_boost && phi_split_slack == other.phi_split_slack;
    }
};

constexpr std::size_t kProfileDepthBuckets = 8;
constexpr std::size_t kProfileSizeBuckets = 8;
constexpr std::size_t kProfilePhiABuckets = 7;
constexpr std::size_t kProfileStageBuckets = 8;
constexpr std::array<u64, kProfileSizeBuckets> kProfileSizeLimits = {
    1'000'000ULL,
    10'000'000ULL,
    100'000'000ULL,
    1'000'000'000ULL,
    10'000'000'000ULL,
    100'000'000'000ULL,
    1'000'000'000'000ULL,
    std::numeric_limits<u64>::max(),
};

std::size_t profile_size_bucket(u64 x) {
    for (std::size_t i = 0; i < kProfileSizeLimits.size(); ++i) {
        if (x <= kProfileSizeLimits[i]) {
            return i;
        }
    }
    return kProfileSizeLimits.size() - 1;
}

const char* profile_size_bucket_label(std::size_t index) {
    static constexpr std::array<const char*, kProfileSizeBuckets> labels = {
        "<=1e6",
        "<=1e7",
        "<=1e8",
        "<=1e9",
        "<=1e10",
        "<=1e11",
        "<=1e12",
        ">1e12",
    };
    return labels[index];
}

std::size_t profile_phi_a_bucket(u64 a) {
    if (a == 0) {
        return 0;
    }
    if (a <= kWheelPrimeCount) {
        return 1;
    }
    if (a <= 16) {
        return 2;
    }
    if (a <= 32) {
        return 3;
    }
    if (a <= 64) {
        return 4;
    }
    if (a <= 128) {
        return 5;
    }
    return 6;
}

const char* profile_phi_a_bucket_label(std::size_t index) {
    static constexpr std::array<const char*, kProfilePhiABuckets> labels = {
        "a=0",
        "a<=6",
        "a<=16",
        "a<=32",
        "a<=64",
        "a<=128",
        "a>128",
    };
    return labels[index];
}

enum class ProfileStage : std::uint8_t {
    Other = 0,
    Prepare,
    TopTerm,
    InnerBound,
    InnerTerm,
    BasePhi,
    SplitPhi,
    RecursiveCount,
};

std::size_t profile_stage_index(ProfileStage stage) {
    return static_cast<std::size_t>(stage);
}

const char* profile_stage_label(std::size_t index) {
    static constexpr std::array<const char*, kProfileStageBuckets> labels = {
        "other",
        "prepare",
        "top_term",
        "inner_bound",
        "inner_term",
        "base_phi",
        "split_phi",
        "recursive_count",
    };
    return labels[index];
}

struct ProfileThreadCounters {
    u64 count_calls = 0;
    u64 count_value_direct = 0;
    u64 count_value_rank_direct = 0;
    u64 count_value_recurse = 0;
    u64 phi_calls = 0;
    u64 phi_cache_hits = 0;
    u64 phi_wheel_hits = 0;
    u64 phi_small_count_hits = 0;
    u64 phi_split_calls = 0;
    u64 phi_split_terms = 0;
    u64 phi_recurrence_calls = 0;
    u32 count_depth = 0;
    ProfileStage stage = ProfileStage::Other;
    std::array<u64, kProfileDepthBuckets> count_calls_by_depth{};
    std::array<u64, kProfileSizeBuckets> count_calls_by_size{};
    std::array<u64, kProfileStageBuckets> count_calls_by_stage{};
    std::array<u64, kProfileStageBuckets> count_value_direct_by_stage{};
    std::array<u64, kProfileStageBuckets> count_value_rank_direct_by_stage{};
    std::array<u64, kProfileStageBuckets> count_value_recurse_by_stage{};
    std::array<u64, kProfileStageBuckets> count_value_outer_recurse_ns_by_stage{};
    std::array<u64, kProfileStageBuckets> phi_calls_by_stage{};
    std::array<u64, kProfilePhiABuckets> phi_calls_by_a{};
    std::array<u64, kProfileSizeBuckets> phi_calls_by_size{};

    void reset() {
        count_calls = 0;
        count_value_direct = 0;
        count_value_rank_direct = 0;
        count_value_recurse = 0;
        phi_calls = 0;
        phi_cache_hits = 0;
        phi_wheel_hits = 0;
        phi_small_count_hits = 0;
        phi_split_calls = 0;
        phi_split_terms = 0;
        phi_recurrence_calls = 0;
        count_depth = 0;
        stage = ProfileStage::Other;
        count_calls_by_depth.fill(0);
        count_calls_by_size.fill(0);
        count_calls_by_stage.fill(0);
        count_value_direct_by_stage.fill(0);
        count_value_rank_direct_by_stage.fill(0);
        count_value_recurse_by_stage.fill(0);
        count_value_outer_recurse_ns_by_stage.fill(0);
        phi_calls_by_stage.fill(0);
        phi_calls_by_a.fill(0);
        phi_calls_by_size.fill(0);
    }
};

thread_local ProfileThreadCounters g_profile_thread;

struct ProfileState {
    bool enabled = false;
    bool detail = false;
    u64 n = 0;
    u64 estimate = 0;
    u64 count_at_estimate = 0;
    u64 result = 0;
    u64 delta = 0;
    u64 top_x = 0;
    u64 top_a = 0;
    u64 top_b = 0;
    u64 top_c = 0;
    u64 top_total = 0;
    u64 top_chunk = 0;
    u32 top_threads = 0;
    u32 config_sieve_limit = 0;
    u32 config_rank_limit = 0;
    u32 config_phi_x = 0;
    u32 config_phi_a = 0;
    u32 config_threads = 0;
    u32 config_a_boost = 0;
    u32 config_phi_split_slack = 0;
    const char* correction_direction = "none";
    double estimate_ms = 0.0;
    double counter_sieve_ms = 0.0;
    double counter_sieve_alloc_ms = 0.0;
    double counter_sieve_mark_ms = 0.0;
    double counter_sieve_pi_ms = 0.0;
    double counter_sieve_rank_ms = 0.0;
    double counter_wheel_ms = 0.0;
    double counter_phi_ms = 0.0;
    double count_setup_ms = 0.0;
    double count_core_ms = 0.0;
    double count_prepare_ms = 0.0;
    double count_base_ms = 0.0;
    double count_parallel_ms = 0.0;
    double count_thread_launch_ms = 0.0;
    double count_join_ms = 0.0;
    double count_worker_min_ms = 0.0;
    double count_worker_avg_ms = 0.0;
    double count_worker_max_ms = 0.0;
    double count_worker_top_sum_ms = 0.0;
    double count_worker_top_avg_ms = 0.0;
    double count_worker_top_max_ms = 0.0;
    double count_worker_phi_split_sum_ms = 0.0;
    double count_worker_phi_split_avg_ms = 0.0;
    double count_worker_phi_split_max_ms = 0.0;
    double count_main_phi_split_ms = 0.0;
    double count_main_top_ms = 0.0;
    double count_partial_reduce_ms = 0.0;
    double correction_ms = 0.0;
    double correction_plan_ms = 0.0;
    double correction_mark_ms = 0.0;
    double correction_scan_ms = 0.0;
    double total_ms = 0.0;
    std::atomic<u64> count_calls{0};
    std::atomic<u64> count_value_direct{0};
    std::atomic<u64> count_value_rank_direct{0};
    std::atomic<u64> count_value_recurse{0};
    std::atomic<u64> phi_calls{0};
    std::atomic<u64> phi_cache_hits{0};
    std::atomic<u64> phi_wheel_hits{0};
    std::atomic<u64> phi_small_count_hits{0};
    std::atomic<u64> phi_split_calls{0};
    std::atomic<u64> phi_split_terms{0};
    std::atomic<u64> phi_recurrence_calls{0};
    std::array<std::atomic<u64>, kProfileDepthBuckets> count_calls_by_depth{};
    std::array<std::atomic<u64>, kProfileSizeBuckets> count_calls_by_size{};
    std::array<std::atomic<u64>, kProfileStageBuckets> count_calls_by_stage{};
    std::array<std::atomic<u64>, kProfileStageBuckets> count_value_direct_by_stage{};
    std::array<std::atomic<u64>, kProfileStageBuckets> count_value_rank_direct_by_stage{};
    std::array<std::atomic<u64>, kProfileStageBuckets> count_value_recurse_by_stage{};
    std::array<std::atomic<u64>, kProfileStageBuckets> count_value_outer_recurse_ns_by_stage{};
    std::array<std::atomic<u64>, kProfileStageBuckets> phi_calls_by_stage{};
    std::array<std::atomic<u64>, kProfilePhiABuckets> phi_calls_by_a{};
    std::array<std::atomic<u64>, kProfileSizeBuckets> phi_calls_by_size{};
    std::atomic<u64> top_iterations{0};
    std::atomic<u64> inner_iterations{0};
    std::atomic<u64> top_phi_split_terms{0};
    std::atomic<u64> correction_chunks{0};
    std::atomic<u64> correction_odds{0};
    std::atomic<u64> correction_mark_entries{0};
    std::atomic<u64> correction_mark_ops{0};

    void reset(bool on, bool detailed = false) {
        enabled = on;
        detail = detailed;
        n = 0;
        estimate = 0;
        count_at_estimate = 0;
        result = 0;
        delta = 0;
        top_x = 0;
        top_a = 0;
        top_b = 0;
        top_c = 0;
        top_total = 0;
        top_chunk = 0;
        top_threads = 0;
        config_sieve_limit = 0;
        config_rank_limit = 0;
        config_phi_x = 0;
        config_phi_a = 0;
        config_threads = 0;
        config_a_boost = 0;
        config_phi_split_slack = 0;
        correction_direction = "none";
        estimate_ms = 0.0;
        counter_sieve_ms = 0.0;
        counter_sieve_alloc_ms = 0.0;
        counter_sieve_mark_ms = 0.0;
        counter_sieve_pi_ms = 0.0;
        counter_sieve_rank_ms = 0.0;
        counter_wheel_ms = 0.0;
        counter_phi_ms = 0.0;
        count_setup_ms = 0.0;
        count_core_ms = 0.0;
        count_prepare_ms = 0.0;
        count_base_ms = 0.0;
        count_parallel_ms = 0.0;
        count_thread_launch_ms = 0.0;
        count_join_ms = 0.0;
        count_worker_min_ms = 0.0;
        count_worker_avg_ms = 0.0;
        count_worker_max_ms = 0.0;
        count_worker_top_sum_ms = 0.0;
        count_worker_top_avg_ms = 0.0;
        count_worker_top_max_ms = 0.0;
        count_worker_phi_split_sum_ms = 0.0;
        count_worker_phi_split_avg_ms = 0.0;
        count_worker_phi_split_max_ms = 0.0;
        count_main_phi_split_ms = 0.0;
        count_main_top_ms = 0.0;
        count_partial_reduce_ms = 0.0;
        correction_ms = 0.0;
        correction_plan_ms = 0.0;
        correction_mark_ms = 0.0;
        correction_scan_ms = 0.0;
        total_ms = 0.0;
        count_calls.store(0, std::memory_order_relaxed);
        count_value_direct.store(0, std::memory_order_relaxed);
        count_value_rank_direct.store(0, std::memory_order_relaxed);
        count_value_recurse.store(0, std::memory_order_relaxed);
        phi_calls.store(0, std::memory_order_relaxed);
        phi_cache_hits.store(0, std::memory_order_relaxed);
        phi_wheel_hits.store(0, std::memory_order_relaxed);
        phi_small_count_hits.store(0, std::memory_order_relaxed);
        phi_split_calls.store(0, std::memory_order_relaxed);
        phi_split_terms.store(0, std::memory_order_relaxed);
        phi_recurrence_calls.store(0, std::memory_order_relaxed);
        for (std::atomic<u64>& value : count_calls_by_depth) {
            value.store(0, std::memory_order_relaxed);
        }
        for (std::atomic<u64>& value : count_calls_by_size) {
            value.store(0, std::memory_order_relaxed);
        }
        for (std::atomic<u64>& value : count_calls_by_stage) {
            value.store(0, std::memory_order_relaxed);
        }
        for (std::atomic<u64>& value : count_value_direct_by_stage) {
            value.store(0, std::memory_order_relaxed);
        }
        for (std::atomic<u64>& value : count_value_rank_direct_by_stage) {
            value.store(0, std::memory_order_relaxed);
        }
        for (std::atomic<u64>& value : count_value_recurse_by_stage) {
            value.store(0, std::memory_order_relaxed);
        }
        for (std::atomic<u64>& value : count_value_outer_recurse_ns_by_stage) {
            value.store(0, std::memory_order_relaxed);
        }
        for (std::atomic<u64>& value : phi_calls_by_stage) {
            value.store(0, std::memory_order_relaxed);
        }
        for (std::atomic<u64>& value : phi_calls_by_a) {
            value.store(0, std::memory_order_relaxed);
        }
        for (std::atomic<u64>& value : phi_calls_by_size) {
            value.store(0, std::memory_order_relaxed);
        }
        top_iterations.store(0, std::memory_order_relaxed);
        inner_iterations.store(0, std::memory_order_relaxed);
        top_phi_split_terms.store(0, std::memory_order_relaxed);
        correction_chunks.store(0, std::memory_order_relaxed);
        correction_odds.store(0, std::memory_order_relaxed);
        correction_mark_entries.store(0, std::memory_order_relaxed);
        correction_mark_ops.store(0, std::memory_order_relaxed);
    }
};

ProfileState g_profile;

void aggregate_profile_thread() {
    if (!g_profile.enabled) {
        return;
    }
    g_profile.count_calls.fetch_add(g_profile_thread.count_calls, std::memory_order_relaxed);
    g_profile.count_value_direct.fetch_add(g_profile_thread.count_value_direct, std::memory_order_relaxed);
    g_profile.count_value_rank_direct.fetch_add(
        g_profile_thread.count_value_rank_direct, std::memory_order_relaxed);
    g_profile.count_value_recurse.fetch_add(g_profile_thread.count_value_recurse, std::memory_order_relaxed);
    g_profile.phi_calls.fetch_add(g_profile_thread.phi_calls, std::memory_order_relaxed);
    g_profile.phi_cache_hits.fetch_add(g_profile_thread.phi_cache_hits, std::memory_order_relaxed);
    g_profile.phi_wheel_hits.fetch_add(g_profile_thread.phi_wheel_hits, std::memory_order_relaxed);
    g_profile.phi_small_count_hits.fetch_add(g_profile_thread.phi_small_count_hits, std::memory_order_relaxed);
    g_profile.phi_split_calls.fetch_add(g_profile_thread.phi_split_calls, std::memory_order_relaxed);
    g_profile.phi_split_terms.fetch_add(g_profile_thread.phi_split_terms, std::memory_order_relaxed);
    g_profile.phi_recurrence_calls.fetch_add(g_profile_thread.phi_recurrence_calls, std::memory_order_relaxed);
    for (std::size_t i = 0; i < kProfileDepthBuckets; ++i) {
        g_profile.count_calls_by_depth[i].fetch_add(
            g_profile_thread.count_calls_by_depth[i], std::memory_order_relaxed);
    }
    for (std::size_t i = 0; i < kProfileSizeBuckets; ++i) {
        g_profile.count_calls_by_size[i].fetch_add(
            g_profile_thread.count_calls_by_size[i], std::memory_order_relaxed);
    }
    for (std::size_t i = 0; i < kProfileStageBuckets; ++i) {
        g_profile.count_calls_by_stage[i].fetch_add(
            g_profile_thread.count_calls_by_stage[i], std::memory_order_relaxed);
        g_profile.count_value_direct_by_stage[i].fetch_add(
            g_profile_thread.count_value_direct_by_stage[i], std::memory_order_relaxed);
        g_profile.count_value_rank_direct_by_stage[i].fetch_add(
            g_profile_thread.count_value_rank_direct_by_stage[i], std::memory_order_relaxed);
        g_profile.count_value_recurse_by_stage[i].fetch_add(
            g_profile_thread.count_value_recurse_by_stage[i], std::memory_order_relaxed);
        g_profile.count_value_outer_recurse_ns_by_stage[i].fetch_add(
            g_profile_thread.count_value_outer_recurse_ns_by_stage[i], std::memory_order_relaxed);
        g_profile.phi_calls_by_stage[i].fetch_add(
            g_profile_thread.phi_calls_by_stage[i], std::memory_order_relaxed);
    }
    for (std::size_t i = 0; i < kProfilePhiABuckets; ++i) {
        g_profile.phi_calls_by_a[i].fetch_add(
            g_profile_thread.phi_calls_by_a[i], std::memory_order_relaxed);
    }
    for (std::size_t i = 0; i < kProfileSizeBuckets; ++i) {
        g_profile.phi_calls_by_size[i].fetch_add(
            g_profile_thread.phi_calls_by_size[i], std::memory_order_relaxed);
    }
    g_profile_thread.reset();
}

struct ProfileCountDepthScope {
    ProfileCountDepthScope() {
        ++g_profile_thread.count_depth;
    }

    ~ProfileCountDepthScope() {
        --g_profile_thread.count_depth;
    }
};

template <bool Profile>
struct ProfileStageScope;

template <>
struct ProfileStageScope<true> {
    ProfileStage previous = ProfileStage::Other;

    explicit ProfileStageScope(ProfileStage next) {
        previous = g_profile_thread.stage;
        g_profile_thread.stage = next;
    }

    ~ProfileStageScope() {
        g_profile_thread.stage = previous;
    }
};

template <>
struct ProfileStageScope<false> {
    explicit ProfileStageScope(ProfileStage next) {
        static_cast<void>(next);
    }
};

CountConfig count_config_for_n(u64 n) {
    if (n >= 24'000'000'000ULL) {
        return {800'000, 800'000, 65'535, 128, 32, 8, 8};
    }
    if (n >= 12'000'000'000ULL) {
        return {500'000, 500'000, 65'535, 128, 24, 8, 8};
    }
    if (n >= 10'000'000'000ULL) {
        return {502'100, 502'100, 131'072, 32, 12, 12, 8};
    }
    if (n >= 9'750'000'000ULL) {
        return {502'100, 1'000'000, 131'072, 32, 12, 12, 8};
    }
    if (n >= 8'000'000'000ULL) {
        return {502'100, 502'100, 131'072, 32, 12, 12, 8};
    }
    if (n >= 4'000'000'000ULL) {
        return {500'000, 500'000, 98'304, 32, 12, 12, 8};
    }
    if (n >= 2'000'000'000ULL) {
        return {260'000, 260'000, 65'535, 32, 12, 12, 8};
    }
    if (n >= 1'500'000'000ULL) {
        return {260'000, 260'000, 36'000, 32, 10, 4, 8};
    }
    if (n >= 800'000'000ULL) {
        return {150'000, 150'000, 28'000, 48, 10, 3, 8};
    }
    if (n >= 200'000'000ULL) {
        return {250'000, 250'000, 12'000, 72, 8, 4, 8};
    }
    if (n >= 50'000'000ULL) {
        return {100'000, 100'000, 2'048, 32, 8, 4, 8};
    }
    if (n >= 20'000'000ULL) {
        return {80'000, 80'000, 2'048, 32, 1, 4, 8};
    }
    if (n < 5'000'000ULL) {
        return {16'000, 16'000, 1'024, 24, 1, 4, 8};
    }
    if (n < 50'000'000ULL) {
        return {32'000, 32'000, 2'048, 32, 1, 4, 8};
    }
    return {750'000, 750'000, 12'000, 72, 1, 4, 8};
}

class LehmerCounter {
public:
    explicit LehmerCounter(CountConfig config) : config_(config) {
        if (g_profile.enabled) {
            const auto sieve_start = std::chrono::steady_clock::now();
            SieveBuildTimings sieve_timings;
            table_.init(config_.sieve_limit, config_.rank_limit, &sieve_timings);
            const auto sieve_stop = std::chrono::steady_clock::now();

            const auto wheel_start = std::chrono::steady_clock::now();
            init_wheel_phi();
            const auto wheel_stop = std::chrono::steady_clock::now();

            const auto phi_start = std::chrono::steady_clock::now();
            init_phi_cache();
            const auto phi_stop = std::chrono::steady_clock::now();

            g_profile.counter_sieve_ms +=
                std::chrono::duration<double, std::milli>(sieve_stop - sieve_start).count();
            g_profile.counter_sieve_alloc_ms += sieve_timings.alloc_ms;
            g_profile.counter_sieve_mark_ms += sieve_timings.mark_ms;
            g_profile.counter_sieve_pi_ms += sieve_timings.pi_ms;
            g_profile.counter_sieve_rank_ms += sieve_timings.rank_ms;
            g_profile.counter_wheel_ms +=
                std::chrono::duration<double, std::milli>(wheel_stop - wheel_start).count();
            g_profile.counter_phi_ms +=
                std::chrono::duration<double, std::milli>(phi_stop - phi_start).count();
        } else {
            table_.init(config_.sieve_limit, config_.rank_limit);
            init_wheel_phi();
            init_phi_cache();
        }
    }

    u64 count(u64 x) {
        return g_profile.enabled ? count_impl<true>(x) : count_impl<false>(x);
    }

    const std::vector<u32>& primes() const {
        return table_.primes;
    }

    u32 limit() const {
        return table_.limit;
    }

    u64 count_parallel(u64 x, u32 thread_count) {
        return g_profile.enabled ?
            count_parallel_impl<true>(x, thread_count) :
            count_parallel_impl<false>(x, thread_count);
    }

    CountConfig config() const {
        return config_;
    }

private:
    CountConfig config_;
    SieveTable table_;
    std::array<u16, 30'030> wheel_phi_prefix_{};
    std::unique_ptr<u16[]> phi_cache_;

    template <bool Profile>
    u64 count_impl(u64 x) {
        if constexpr (Profile) {
            ++g_profile_thread.count_calls;
            ++g_profile_thread.count_calls_by_depth[
                std::min<std::size_t>(g_profile_thread.count_depth, kProfileDepthBuckets - 1)];
            ++g_profile_thread.count_calls_by_size[profile_size_bucket(x)];
            ++g_profile_thread.count_calls_by_stage[profile_stage_index(g_profile_thread.stage)];
            ProfileCountDepthScope depth_scope;
            return count_impl_body<true>(x);
        }
        return count_impl_body<false>(x);
    }

    template <bool Profile>
    u64 count_impl_body(u64 x) {
        if (x <= table_.limit) {
            return table_.count(x);
        }

        u64 base_a = 0;
        u64 b = 0;
        u64 c = 0;
        {
            ProfileStageScope<Profile> stage(ProfileStage::Prepare);
            base_a = count_value<Profile>(iroot4(x));
            b = count_value<Profile>(isqrt(x));
            c = count_value<Profile>(icbrt(x));
        }
        const u64 a = std::min<u64>(c, base_a + config_.a_boost);

        u64 sum = 0;
        {
            ProfileStageScope<Profile> stage(ProfileStage::BasePhi);
            sum = phi<Profile>(x, a) + ((b + a - 2) * (b - a + 1)) / 2;
        }
        for (u64 i = a; i < b; ++i) {
            const u64 w = x / table_.primes[static_cast<std::size_t>(i)];
            {
                ProfileStageScope<Profile> stage(ProfileStage::TopTerm);
                sum -= count_value<Profile>(w);
            }

            if (i < c) {
                u64 bi = 0;
                {
                    ProfileStageScope<Profile> stage(ProfileStage::InnerBound);
                    bi = count_value<Profile>(isqrt(w));
                }
                for (u64 j = i; j < bi; ++j) {
                    ProfileStageScope<Profile> stage(ProfileStage::InnerTerm);
                    sum -= count_value<Profile>(w / table_.primes[static_cast<std::size_t>(j)]) - j;
                }
            }
        }

        return sum;
    }

    template <bool Profile>
    u64 count_parallel_impl(u64 x, u32 thread_count) {
        if (x <= table_.limit || thread_count <= 1) {
            return count_impl<Profile>(x);
        }

        std::chrono::steady_clock::time_point prepare_start;
        if constexpr (Profile) {
            prepare_start = std::chrono::steady_clock::now();
        }
        u64 base_a = 0;
        u64 b = 0;
        u64 c = 0;
        {
            ProfileStageScope<Profile> stage(ProfileStage::Prepare);
            base_a = count_value<Profile>(iroot4(x));
            b = count_value<Profile>(isqrt(x));
            c = count_value<Profile>(icbrt(x));
        }
        const u64 a = std::min<u64>(c, base_a + config_.a_boost);
        const u64 arithmetic = ((b + a - 2) * (b - a + 1)) / 2;
        const u64 total = b - a;
        const bool split_base_phi = a > config_.phi_a + config_.phi_split_slack;
        const u64 base_phi_a = split_base_phi && config_.phi_a > kWheelPrimeCount ?
            config_.phi_a - 1 :
            (split_base_phi ? kWheelPrimeCount : a);
        if (total < 4'096) {
            if constexpr (Profile) {
                g_profile.count_prepare_ms +=
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - prepare_start)
                        .count();
            }
            return count_impl<Profile>(x);
        }
        if constexpr (Profile) {
            g_profile.count_prepare_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - prepare_start)
                    .count();
        }

        thread_count = static_cast<u32>(std::min<u64>(thread_count, total));
        std::vector<u64> partials(thread_count, 0);
        std::vector<double> worker_ms;
        std::vector<double> worker_top_ms;
        std::vector<double> worker_phi_split_ms;
        if constexpr (Profile) {
            worker_ms.assign(thread_count, 0.0);
            worker_top_ms.assign(thread_count, 0.0);
            worker_phi_split_ms.assign(thread_count, 0.0);
        }
        std::vector<std::thread> threads;
        threads.reserve(thread_count);

        std::atomic<u64> next_i{a};
        std::atomic<u64> next_phi_i{base_phi_a};
        const u64 chunk = total > 30'000 ?
            NPRIME_TOP_CHUNK_BIG :
            (total > 20'000 ? NPRIME_TOP_CHUNK_MID : NPRIME_TOP_CHUNK_SMALL);
        if constexpr (Profile) {
            g_profile.top_x = x;
            g_profile.top_a = a;
            g_profile.top_b = b;
            g_profile.top_c = c;
            g_profile.top_total = total;
            g_profile.top_chunk = chunk;
            g_profile.top_threads = thread_count;
            if (split_base_phi) {
                g_profile.phi_split_calls.fetch_add(1, std::memory_order_relaxed);
                g_profile.phi_split_terms.fetch_add(a - base_phi_a, std::memory_order_relaxed);
                g_profile.top_phi_split_terms.fetch_add(a - base_phi_a, std::memory_order_relaxed);
            }
        }
        std::chrono::steady_clock::time_point parallel_start;
        if constexpr (Profile) {
            parallel_start = std::chrono::steady_clock::now();
        }
        std::chrono::steady_clock::time_point launch_start;
        if constexpr (Profile) {
            launch_start = std::chrono::steady_clock::now();
        }
        for (u32 t = 0; t < thread_count; ++t) {
            threads.emplace_back([this, x, a, c, b, split_base_phi, chunk,
                                  &next_i, &next_phi_i, &partials, &worker_ms,
                                  &worker_top_ms, &worker_phi_split_ms, t]() {
                static_cast<void>(worker_ms);
                static_cast<void>(worker_top_ms);
                static_cast<void>(worker_phi_split_ms);
                if constexpr (Profile) {
                    g_profile_thread.reset();
                }
                std::chrono::steady_clock::time_point worker_start;
                std::chrono::steady_clock::time_point stage_start;
                if constexpr (Profile) {
                    worker_start = std::chrono::steady_clock::now();
                    stage_start = worker_start;
                }
                u64 local = 0;
                u64 local_top_iterations = 0;
                u64 local_inner_iterations = 0;
                for (;;) {
                    const u64 begin = next_i.fetch_add(chunk, std::memory_order_relaxed);
                    if (begin >= b) {
                        break;
                    }
                    const u64 end = std::min<u64>(b, begin + chunk);
                    if constexpr (Profile) {
                        local_top_iterations += end - begin;
                    }
                    for (u64 i = begin; i < end; ++i) {
                        const u64 w = x / table_.primes[static_cast<std::size_t>(i)];
                        {
                            ProfileStageScope<Profile> stage(ProfileStage::TopTerm);
                            local += count_value<Profile>(w);
                        }

                        if (i < c) {
                            u64 bi = 0;
                            {
                                ProfileStageScope<Profile> stage(ProfileStage::InnerBound);
                                bi = count_value<Profile>(isqrt(w));
                            }
                            for (u64 j = i; j < bi; ++j) {
                                if constexpr (Profile) {
                                    ++local_inner_iterations;
                                }
                                ProfileStageScope<Profile> stage(ProfileStage::InnerTerm);
                                local += count_value<Profile>(w / table_.primes[static_cast<std::size_t>(j)]) - j;
                            }
                        }
                    }
                }
                if constexpr (Profile) {
                    const auto top_stop = std::chrono::steady_clock::now();
                    worker_top_ms[t] =
                        std::chrono::duration<double, std::milli>(top_stop - stage_start).count();
                    stage_start = top_stop;
                }
                if (split_base_phi) {
                    ProfileStageScope<Profile> stage(ProfileStage::SplitPhi);
                    for (;;) {
                        const u64 begin = next_phi_i.fetch_add(chunk, std::memory_order_relaxed);
                        if (begin >= a) {
                            break;
                        }
                        const u64 end = std::min<u64>(a, begin + chunk);
                        for (u64 i = begin; i < end; ++i) {
                            local += phi<Profile>(x / table_.primes[static_cast<std::size_t>(i)], i);
                        }
                    }
                }
                if constexpr (Profile) {
                    if (split_base_phi) {
                        worker_phi_split_ms[t] =
                            std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - stage_start)
                                .count();
                    }
                    g_profile.top_iterations.fetch_add(local_top_iterations, std::memory_order_relaxed);
                    g_profile.inner_iterations.fetch_add(local_inner_iterations, std::memory_order_relaxed);
                    aggregate_profile_thread();
                    worker_ms[t] =
                        std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - worker_start)
                            .count();
                }
                partials[t] = local;
            });
        }
        if constexpr (Profile) {
            g_profile.count_thread_launch_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - launch_start)
                    .count();
        }

        std::chrono::steady_clock::time_point base_start;
        if constexpr (Profile) {
            base_start = std::chrono::steady_clock::now();
        }
        u64 base = 0;
        {
            ProfileStageScope<Profile> stage(ProfileStage::BasePhi);
            base = phi<Profile>(x, base_phi_a) + arithmetic;
        }
        if constexpr (Profile) {
            g_profile.count_base_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - base_start)
                    .count();
        }

        u64 main_subtract = 0;
        std::chrono::steady_clock::time_point main_phi_split_start;
        if constexpr (Profile) {
            main_phi_split_start = std::chrono::steady_clock::now();
        }
        if (split_base_phi) {
            ProfileStageScope<Profile> stage(ProfileStage::SplitPhi);
            for (;;) {
                const u64 begin = next_phi_i.fetch_add(chunk, std::memory_order_relaxed);
                if (begin >= a) {
                    break;
                }
                const u64 end = std::min<u64>(a, begin + chunk);
                for (u64 i = begin; i < end; ++i) {
                    main_subtract += phi<Profile>(x / table_.primes[static_cast<std::size_t>(i)], i);
                }
            }
        }
        if constexpr (Profile) {
            if (split_base_phi) {
                g_profile.count_main_phi_split_ms +=
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - main_phi_split_start)
                        .count();
            }
        }
        u64 main_top_iterations = 0;
        u64 main_inner_iterations = 0;
        std::chrono::steady_clock::time_point main_top_start;
        if constexpr (Profile) {
            main_top_start = std::chrono::steady_clock::now();
        }
        for (;;) {
            const u64 begin = next_i.fetch_add(chunk, std::memory_order_relaxed);
            if (begin >= b) {
                break;
            }
            const u64 end = std::min<u64>(b, begin + chunk);
            if constexpr (Profile) {
                main_top_iterations += end - begin;
            }
            for (u64 i = begin; i < end; ++i) {
                const u64 w = x / table_.primes[static_cast<std::size_t>(i)];
                {
                    ProfileStageScope<Profile> stage(ProfileStage::TopTerm);
                    main_subtract += count_value<Profile>(w);
                }

                if (i < c) {
                    u64 bi = 0;
                    {
                        ProfileStageScope<Profile> stage(ProfileStage::InnerBound);
                        bi = count_value<Profile>(isqrt(w));
                    }
                    for (u64 j = i; j < bi; ++j) {
                        if constexpr (Profile) {
                            ++main_inner_iterations;
                        }
                        ProfileStageScope<Profile> stage(ProfileStage::InnerTerm);
                        main_subtract += count_value<Profile>(w / table_.primes[static_cast<std::size_t>(j)]) - j;
                    }
                }
            }
        }
        if constexpr (Profile) {
            g_profile.count_main_top_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - main_top_start)
                    .count();
            g_profile.top_iterations.fetch_add(main_top_iterations, std::memory_order_relaxed);
            g_profile.inner_iterations.fetch_add(main_inner_iterations, std::memory_order_relaxed);
        }

        u64 subtract = main_subtract;
        std::chrono::steady_clock::time_point join_start;
        if constexpr (Profile) {
            join_start = std::chrono::steady_clock::now();
        }
        for (std::thread& thread : threads) {
            thread.join();
        }
        if constexpr (Profile) {
            g_profile.count_join_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - join_start)
                    .count();
            double worker_sum = 0.0;
            double worker_min = std::numeric_limits<double>::max();
            double worker_max = 0.0;
            double worker_top_sum = 0.0;
            double worker_top_max = 0.0;
            double worker_phi_split_sum = 0.0;
            double worker_phi_split_max = 0.0;
            for (const double ms : worker_ms) {
                worker_sum += ms;
                worker_min = std::min(worker_min, ms);
                worker_max = std::max(worker_max, ms);
            }
            for (const double ms : worker_top_ms) {
                worker_top_sum += ms;
                worker_top_max = std::max(worker_top_max, ms);
            }
            for (const double ms : worker_phi_split_ms) {
                worker_phi_split_sum += ms;
                worker_phi_split_max = std::max(worker_phi_split_max, ms);
            }
            g_profile.count_worker_min_ms = worker_min;
            g_profile.count_worker_avg_ms = worker_sum / static_cast<double>(worker_ms.size());
            g_profile.count_worker_max_ms = worker_max;
            g_profile.count_worker_top_sum_ms = worker_top_sum;
            g_profile.count_worker_top_avg_ms = worker_top_sum / static_cast<double>(worker_top_ms.size());
            g_profile.count_worker_top_max_ms = worker_top_max;
            g_profile.count_worker_phi_split_sum_ms = worker_phi_split_sum;
            g_profile.count_worker_phi_split_avg_ms =
                worker_phi_split_sum / static_cast<double>(worker_phi_split_ms.size());
            g_profile.count_worker_phi_split_max_ms = worker_phi_split_max;
        }
        std::chrono::steady_clock::time_point reduce_start;
        if constexpr (Profile) {
            reduce_start = std::chrono::steady_clock::now();
        }
        for (u64 part : partials) {
            subtract += part;
        }
        if constexpr (Profile) {
            g_profile.count_partial_reduce_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - reduce_start)
                    .count();
        }
        if constexpr (Profile) {
            g_profile.count_parallel_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - parallel_start)
                    .count();
        }

        return base - subtract;
    }

    template <bool Profile>
    u64 count_value(u64 x) {
        if (x <= table_.limit) {
            if constexpr (Profile) {
                ++g_profile_thread.count_value_direct;
                ++g_profile_thread.count_value_direct_by_stage[profile_stage_index(g_profile_thread.stage)];
            }
            return table_.count(x);
        }
        if (x <= table_.rank_limit) {
            if constexpr (Profile) {
                ++g_profile_thread.count_value_rank_direct;
                ++g_profile_thread.count_value_direct;
                ++g_profile_thread.count_value_rank_direct_by_stage[profile_stage_index(g_profile_thread.stage)];
                ++g_profile_thread.count_value_direct_by_stage[profile_stage_index(g_profile_thread.stage)];
            }
            return table_.rank_count(x);
        }
        if constexpr (Profile) {
            const std::size_t caller_stage = profile_stage_index(g_profile_thread.stage);
            const bool outer_recursive_call = g_profile_thread.count_depth <= 1;
            std::chrono::steady_clock::time_point recurse_start;
            if (outer_recursive_call) {
                recurse_start = std::chrono::steady_clock::now();
            }
            ++g_profile_thread.count_value_recurse;
            ++g_profile_thread.count_value_recurse_by_stage[caller_stage];
            {
                ProfileStageScope<true> stage(ProfileStage::RecursiveCount);
                const u64 result = count_impl<true>(x);
                if (outer_recursive_call) {
                    g_profile_thread.count_value_outer_recurse_ns_by_stage[caller_stage] +=
                        elapsed_ns(recurse_start, std::chrono::steady_clock::now());
                }
                return result;
            }
        }
        return count_impl<false>(x);
    }

    void init_wheel_phi() {
        u16 count = 0;
        for (u32 x = 0; x < wheel_phi_prefix_.size(); ++x) {
            if (x > 0 && x % 2 != 0 && x % 3 != 0 && x % 5 != 0 &&
                x % 7 != 0 && x % 11 != 0 && x % 13 != 0) {
                ++count;
            }
            wheel_phi_prefix_[x] = count;
        }
    }

    u64 wheel_phi(u64 x) const {
        static constexpr u64 kWheel = 30'030;
        static constexpr u64 kWheelCount = 5'760;
        return (x / kWheel) * kWheelCount + wheel_phi_prefix_[static_cast<std::size_t>(x % kWheel)];
    }

    u64 wheel_phi7(u64 x) const {
        return wheel_phi(x) - wheel_phi(x / 17);
    }

    u64 wheel_phi8(u64 x) const {
        return wheel_phi7(x) - wheel_phi7(x / 19);
    }

    u64 wheel_phi9(u64 x) const {
        return wheel_phi8(x) - wheel_phi8(x / 23);
    }

    u64 wheel_phi10(u64 x) const {
        return wheel_phi9(x) - wheel_phi9(x / 29);
    }

    void init_phi_cache() {
        if (config_.phi_a <= kWheelPrimeCount) {
            phi_cache_.reset();
            return;
        }
        const std::size_t rows = static_cast<std::size_t>(config_.phi_a - kWheelPrimeCount);
        phi_cache_.reset(new u16[rows * config_.phi_x]);
        for (u32 x = 0; x < config_.phi_x; ++x) {
            const u64 value = wheel_phi(x);
            if (value > std::numeric_limits<u16>::max()) {
                throw std::overflow_error("phi cache wheel row exceeds uint16 range");
            }
            phi_cache_[x] = static_cast<u16>(value);
        }

        for (u32 a = kWheelPrimeCount + 1; a < config_.phi_a; ++a) {
            const u32 p = table_.primes[static_cast<std::size_t>(a - 1)];
            const std::size_t row = static_cast<std::size_t>(a - kWheelPrimeCount) * config_.phi_x;
            const std::size_t prev = static_cast<std::size_t>(a - kWheelPrimeCount - 1) * config_.phi_x;
            for (u32 quotient = 0, begin = 0; begin < config_.phi_x; ++quotient) {
                const u32 end = std::min<u32>(config_.phi_x, begin + p);
                const u16 remove = phi_cache_[prev + quotient];
                for (u32 x = begin; x < end; ++x) {
                    phi_cache_[row + x] = static_cast<u16>(phi_cache_[prev + x] - remove);
                }
                begin = end;
            }
        }
    }

    template <bool Profile>
    u64 phi(u64 x, u64 a) {
        if constexpr (Profile) {
            ++g_profile_thread.phi_calls;
            ++g_profile_thread.phi_calls_by_stage[profile_stage_index(g_profile_thread.stage)];
            if (g_profile.detail) {
                ++g_profile_thread.phi_calls_by_a[profile_phi_a_bucket(a)];
                ++g_profile_thread.phi_calls_by_size[profile_size_bucket(x)];
            }
        }
        if (a == 0) {
            return x;
        }
        if (a >= kWheelPrimeCount && a < config_.phi_a && x < config_.phi_x) {
            if constexpr (Profile) {
                ++g_profile_thread.phi_cache_hits;
            }
            const std::size_t row = static_cast<std::size_t>(a - kWheelPrimeCount) * config_.phi_x;
            return phi_cache_[row + static_cast<std::size_t>(x)];
        }
        if (a == kWheelPrimeCount) {
            if constexpr (Profile) {
                ++g_profile_thread.phi_wheel_hits;
            }
            return wheel_phi(x);
        }
        if (x <= table_.limit && a < table_.primes.size()) {
            const u64 next = table_.primes[static_cast<std::size_t>(a)];
            if (x < next * next) {
                if constexpr (Profile) {
                    ++g_profile_thread.phi_small_count_hits;
                }
                return table_.count(x) - a + 1;
            }
        }
#if NPRIME_DIRECT_WHEEL_A >= 7
        if (a == 7) {
            if constexpr (Profile) {
                ++g_profile_thread.phi_wheel_hits;
            }
            return wheel_phi7(x);
        }
#endif
#if NPRIME_DIRECT_WHEEL_A >= 8
        if (a == 8) {
            if constexpr (Profile) {
                ++g_profile_thread.phi_wheel_hits;
            }
            return wheel_phi8(x);
        }
#endif
#if NPRIME_DIRECT_WHEEL_A >= 9
        if (a == 9) {
            if constexpr (Profile) {
                ++g_profile_thread.phi_wheel_hits;
            }
            return wheel_phi9(x);
        }
#endif
#if NPRIME_DIRECT_WHEEL_A >= 10
        if (a == 10) {
            if constexpr (Profile) {
                ++g_profile_thread.phi_wheel_hits;
            }
            return wheel_phi10(x);
        }
#endif
        u64 result = 0;
        if (a > config_.phi_a + config_.phi_split_slack) {
            const u64 b = config_.phi_a > kWheelPrimeCount ? config_.phi_a - 1 : kWheelPrimeCount;
            if constexpr (Profile) {
                ++g_profile_thread.phi_split_calls;
                g_profile_thread.phi_split_terms += a - b;
            }
            result = phi<Profile>(x, b);
            for (u64 i = b; i < a; ++i) {
                result -= phi<Profile>(x / table_.primes[static_cast<std::size_t>(i)], i);
            }
        } else {
            if constexpr (Profile) {
                ++g_profile_thread.phi_recurrence_calls;
            }
            result = phi<Profile>(x, a - 1) -
                phi<Profile>(x / table_.primes[static_cast<std::size_t>(a - 1)], a - 1);
        }
        return result;
    }
};

std::unique_ptr<LehmerCounter>& lehmer_counter_storage() {
    static std::unique_ptr<LehmerCounter> counter;
    return counter;
}

LehmerCounter& lehmer_counter(CountConfig config) {
    std::unique_ptr<LehmerCounter>& counter = lehmer_counter_storage();
    if (!counter || !(counter->config() == config)) {
        counter = std::make_unique<LehmerCounter>(config);
    }
    return *counter;
}

u64 nth_prime_estimate(u64 n) {
    switch (n) {
    case 0:
        return 0;
    case 1:
        return 2;
    case 2:
        return 3;
    case 3:
        return 5;
    case 4:
        return 7;
    case 5:
        return 11;
    default:
        break;
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

long double riemann_r(long double x) {
    if (x < 2.0L) {
        return 0.0L;
    }

    static constexpr long double kZeta[] = {
        1.64493406684822643647L,
        1.20205690315959428540L,
        1.08232323371113819152L,
        1.03692775514336992633L,
        1.01734306198444913971L,
        1.00834927738192282684L,
        1.00407735619794433938L,
        1.00200839282608221442L,
        1.00099457512781808534L,
        1.00049418860411946456L,
        1.00024608655330804830L,
        1.00012271334757848915L,
        1.00006124813505870483L,
        1.00003058823630702049L,
        1.00001528225940865187L,
        1.00000763719763789976L,
        1.00000381729326499984L,
        1.00000190821271655394L,
        1.00000095396203387280L,
        1.00000047693298678781L,
        1.00000023845050272773L,
        1.00000011921992596531L,
        1.00000005960818905126L,
        1.00000002980350351465L,
    };

    const long double logx = std::log(x);
    long double power = 1.0L;
    long double factorial = 1.0L;
    long double sum = 1.0L;

    for (u32 k = 1; k <= 60; ++k) {
        power *= logx;
        factorial *= static_cast<long double>(k);
        const long double zeta = k <= (sizeof(kZeta) / sizeof(kZeta[0])) ? kZeta[k - 1] : 1.0L;
        sum += power / (static_cast<long double>(k) * factorial * zeta);
    }

    return sum;
}

u64 nth_prime_riemann_estimate(u64 n) {
    if (n < 100'000) {
        return nth_prime_estimate(n);
    }

    long double x = static_cast<long double>(nth_prime_estimate(n));
    const long double target = static_cast<long double>(n);
    for (u32 i = 0; i < 4; ++i) {
        x -= (riemann_r(x) - target) * std::log(x);
        if (x < 2.0L) {
            x = 2.0L;
        }
    }

    return static_cast<u64>(x);
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

thread_local CountConfig g_count_config =
    {NPRIME_SIEVE_LIMIT, NPRIME_SIEVE_LIMIT, NPRIME_PHI_X, NPRIME_PHI_A, 1, 4, 8};
thread_local bool g_count_config_forced = false;
thread_local CountConfig g_forced_count_config =
    {NPRIME_SIEVE_LIMIT, NPRIME_SIEVE_LIMIT, NPRIME_PHI_X, NPRIME_PHI_A, 1, 4, 8};

struct ForcedCountConfigScope {
    bool previous_forced = false;
    CountConfig previous_config =
        {NPRIME_SIEVE_LIMIT, NPRIME_SIEVE_LIMIT, NPRIME_PHI_X, NPRIME_PHI_A, 1, 4, 8};

    explicit ForcedCountConfigScope(CountConfig config) {
        previous_forced = g_count_config_forced;
        previous_config = g_forced_count_config;
        g_forced_count_config = config;
        g_count_config_forced = true;
    }

    ~ForcedCountConfigScope() {
        g_forced_count_config = previous_config;
        g_count_config_forced = previous_forced;
    }
};

CountConfig selected_count_config(u64 n) {
    return g_count_config_forced ? g_forced_count_config : count_config_for_n(n);
}

void ensure_count_config_covers_estimate(u64 n, u64 estimate) {
    const u64 needed_limit = isqrt(estimate + initial_radius(n)) + 1;
    if (needed_limit > std::numeric_limits<u32>::max()) {
        throw std::overflow_error("required sieve limit exceeds uint32 range");
    }
    if (needed_limit > g_count_config.sieve_limit) {
        g_count_config.sieve_limit = static_cast<u32>(needed_limit);
    }
    if (g_count_config.rank_limit < g_count_config.sieve_limit) {
        g_count_config.rank_limit = g_count_config.sieve_limit;
    }
}

u64 prime_pi_with_timing(u64 x, double* setup_ms, double* core_ms) {
    const auto setup_start = std::chrono::steady_clock::now();
    LehmerCounter& counter = lehmer_counter(g_count_config);
    const auto setup_stop = std::chrono::steady_clock::now();

    u64 result = 0;
    const auto core_start = std::chrono::steady_clock::now();
    if (g_count_config.threads > 1) {
        const u32 threads = std::clamp<u32>(g_count_config.threads, 2, 32);
        result = counter.count_parallel(x, threads);
    } else {
        result = counter.count(x);
    }
    const auto core_stop = std::chrono::steady_clock::now();

    if (setup_ms != nullptr) {
        *setup_ms = std::chrono::duration<double, std::milli>(setup_stop - setup_start).count();
    }
    if (core_ms != nullptr) {
        *core_ms = std::chrono::duration<double, std::milli>(core_stop - core_start).count();
    }
    if (g_profile.enabled) {
        g_profile.count_setup_ms =
            std::chrono::duration<double, std::milli>(setup_stop - setup_start).count();
        g_profile.count_core_ms =
            std::chrono::duration<double, std::milli>(core_stop - core_start).count();
    }
    return result;
}

u64 prime_pi(u64 x) {
    return prime_pi_with_timing(x, nullptr, nullptr);
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

thread_local SegmentPlan g_segment_plan;

const std::vector<SegmentEntry>& segment_plan(const std::vector<u32>& primes, u64 low, u64 high, u64 odds) {
    SegmentPlan& plan = g_segment_plan;
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

struct SearchState {
    u64 n = 0;
    u64 estimate = 0;
    u64 count_at_estimate = 0;
    bool table_primes_ok = false;
    bool valid = false;
};

thread_local SearchState g_cached_search;

void reset_runtime_state() {
    lehmer_counter_storage().reset();
    g_segment_plan = SegmentPlan{};
    g_cached_search = SearchState{};
    g_count_config = {NPRIME_SIEVE_LIMIT, NPRIME_SIEVE_LIMIT, NPRIME_PHI_X, NPRIME_PHI_A, 1, 4, 8};
}

template <bool Profile>
u64 nth_prime_after_impl(const std::vector<u32>& primes, u64 start, u64 offset) {
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
        if constexpr (Profile) {
            g_profile.correction_chunks.fetch_add(1, std::memory_order_relaxed);
            g_profile.correction_odds.fetch_add(chunk_odds, std::memory_order_relaxed);
        }
        std::array<std::uint8_t, kStackOddChunk> stack_composite;
        std::vector<std::uint8_t> composite;
        std::uint8_t* composite_data = stack_composite.data();
        if (chunk_odds <= kStackOddChunk) {
            std::fill_n(composite_data, static_cast<std::size_t>(chunk_odds), 0);
        } else {
            composite.assign(static_cast<std::size_t>(chunk_odds), 0);
            composite_data = composite.data();
        }

        std::chrono::steady_clock::time_point stage_start;
        if constexpr (Profile) {
            stage_start = std::chrono::steady_clock::now();
        }
        const std::vector<SegmentEntry>& entries = segment_plan(primes, chunk_low, chunk_high, chunk_odds);
        if constexpr (Profile) {
            g_profile.correction_plan_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - stage_start)
                    .count();
            stage_start = std::chrono::steady_clock::now();
        }
        for (const SegmentEntry& entry : entries) {
            if constexpr (Profile) {
                g_profile.correction_mark_entries.fetch_add(1, std::memory_order_relaxed);
            }
            for (u32 index = entry.first_index; index < chunk_odds; index += entry.step) {
                if constexpr (Profile) {
                    g_profile.correction_mark_ops.fetch_add(1, std::memory_order_relaxed);
                }
                composite_data[index] = 1;
            }
        }
        if constexpr (Profile) {
            g_profile.correction_mark_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - stage_start)
                    .count();
            stage_start = std::chrono::steady_clock::now();
        }

        for (u64 i = 0; i < chunk_odds; ++i) {
            if (!composite_data[static_cast<std::size_t>(i)]) {
                ++seen;
                if (seen == offset) {
                    if constexpr (Profile) {
                        g_profile.correction_scan_ms +=
                            std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - stage_start)
                                .count();
                    }
                    return chunk_low + 2 * i;
                }
            }
        }
        if constexpr (Profile) {
            g_profile.correction_scan_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - stage_start)
                    .count();
        }

        if (chunk_high > std::numeric_limits<u64>::max() - 2) {
            break;
        }
        chunk_low = chunk_high + 2;
    }

    throw std::overflow_error("nth prime exceeds uint64 range");
}

template <bool Profile>
u64 nth_prime_before_or_at_impl(const std::vector<u32>& primes, u64 start, u64 skip) {
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
        if constexpr (Profile) {
            g_profile.correction_chunks.fetch_add(1, std::memory_order_relaxed);
            g_profile.correction_odds.fetch_add(chunk_odds, std::memory_order_relaxed);
        }
        std::array<std::uint8_t, kStackOddChunk> stack_composite;
        std::vector<std::uint8_t> composite;
        std::uint8_t* composite_data = stack_composite.data();
        if (chunk_odds <= kStackOddChunk) {
            std::fill_n(composite_data, static_cast<std::size_t>(chunk_odds), 0);
        } else {
            composite.assign(static_cast<std::size_t>(chunk_odds), 0);
            composite_data = composite.data();
        }

        std::chrono::steady_clock::time_point stage_start;
        if constexpr (Profile) {
            stage_start = std::chrono::steady_clock::now();
        }
        const std::vector<SegmentEntry>& entries = segment_plan(primes, chunk_low, chunk_high, chunk_odds);
        if constexpr (Profile) {
            g_profile.correction_plan_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - stage_start)
                    .count();
            stage_start = std::chrono::steady_clock::now();
        }
        for (const SegmentEntry& entry : entries) {
            if constexpr (Profile) {
                g_profile.correction_mark_entries.fetch_add(1, std::memory_order_relaxed);
            }
            for (u32 index = entry.first_index; index < chunk_odds; index += entry.step) {
                if constexpr (Profile) {
                    g_profile.correction_mark_ops.fetch_add(1, std::memory_order_relaxed);
                }
                composite_data[index] = 1;
            }
        }
        if constexpr (Profile) {
            g_profile.correction_mark_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - stage_start)
                    .count();
            stage_start = std::chrono::steady_clock::now();
        }

        for (u64 i = chunk_odds; i-- > 0;) {
            if (!composite_data[static_cast<std::size_t>(i)]) {
                if (seen == skip) {
                    if constexpr (Profile) {
                        g_profile.correction_scan_ms +=
                            std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - stage_start)
                                .count();
                    }
                    return chunk_low + 2 * i;
                }
                ++seen;
            }
        }
        if constexpr (Profile) {
            g_profile.correction_scan_ms +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - stage_start)
                    .count();
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
    g_count_config = selected_count_config(n);
    if (g_profile.enabled) {
        g_profile.n = n;
        g_profile.config_sieve_limit = g_count_config.sieve_limit;
        g_profile.config_rank_limit = g_count_config.rank_limit;
        g_profile.config_phi_x = g_count_config.phi_x;
        g_profile.config_phi_a = g_count_config.phi_a;
        g_profile.config_threads = g_count_config.threads;
        g_profile.config_a_boost = g_count_config.a_boost;
        g_profile.config_phi_split_slack = g_count_config.phi_split_slack;
    }

    u64 estimate = 0;
    u64 count_at_estimate = 0;
    bool table_primes_ok = false;
    if (g_cached_search.valid && g_cached_search.n == n) {
        estimate = g_cached_search.estimate;
        ensure_count_config_covers_estimate(n, estimate);
        count_at_estimate = g_cached_search.count_at_estimate;
        table_primes_ok = g_cached_search.table_primes_ok;
    } else {
        const auto estimate_start = std::chrono::steady_clock::now();
        estimate = nth_prime_riemann_estimate(n);
        const auto estimate_stop = std::chrono::steady_clock::now();
        if (g_profile.enabled) {
            g_profile.estimate_ms =
                std::chrono::duration<double, std::milli>(estimate_stop - estimate_start).count();
        }
        ensure_count_config_covers_estimate(n, estimate);
        if (g_profile.enabled) {
            g_profile.config_sieve_limit = g_count_config.sieve_limit;
            g_profile.config_rank_limit = g_count_config.rank_limit;
        }
        count_at_estimate = prime_pi(estimate);
        if (g_profile.enabled) {
            g_profile.estimate = estimate;
            g_profile.count_at_estimate = count_at_estimate;
            g_profile.delta = count_at_estimate > n ? count_at_estimate - n : n - count_at_estimate;
        }

        const u64 needed_limit = isqrt(estimate + initial_radius(n)) + 1;
        table_primes_ok = needed_limit <= lehmer_counter(g_count_config).limit();
        g_cached_search = {n, estimate, count_at_estimate, table_primes_ok, true};
    }

    const std::vector<u32>* primes = &lehmer_counter(g_count_config).primes();
    std::vector<u32> generated_primes;
    if (!table_primes_ok) {
        const u64 needed_limit = isqrt(estimate + initial_radius(n)) + 1;
        generated_primes = simple_primes(needed_limit);
        primes = &generated_primes;
    }

    if (count_at_estimate < n) {
        if (g_profile.enabled) {
            g_profile.correction_direction = "forward";
        }
        const auto correction_start = std::chrono::steady_clock::now();
        const u64 result = g_profile.enabled ?
            nth_prime_after_impl<true>(*primes, estimate, n - count_at_estimate) :
            nth_prime_after_impl<false>(*primes, estimate, n - count_at_estimate);
        const auto correction_stop = std::chrono::steady_clock::now();
        if (g_profile.enabled) {
            g_profile.correction_ms =
                std::chrono::duration<double, std::milli>(correction_stop - correction_start).count();
        }
        return result;
    }
    if (g_profile.enabled) {
        g_profile.correction_direction = "backward";
    }
    const auto correction_start = std::chrono::steady_clock::now();
    const u64 result = g_profile.enabled ?
        nth_prime_before_or_at_impl<true>(*primes, estimate, count_at_estimate - n) :
        nth_prime_before_or_at_impl<false>(*primes, estimate, count_at_estimate - n);
    const auto correction_stop = std::chrono::steady_clock::now();
    if (g_profile.enabled) {
        g_profile.correction_ms =
            std::chrono::duration<double, std::milli>(correction_stop - correction_start).count();
    }
    return result;
}

void trace_nth_prime(u64 n) {
    if (n == 0) {
        std::cout << "n=0 result=0\n";
        return;
    }

    g_count_config = selected_count_config(n);
    u64 estimate = nth_prime_riemann_estimate(n);
    ensure_count_config_covers_estimate(n, estimate);
    double setup_ms = 0.0;
    double core_ms = 0.0;
    u64 count_at_estimate = prime_pi_with_timing(estimate, &setup_ms, &core_ms);
    const u64 initial_delta = count_at_estimate > n ? count_at_estimate - n : n - count_at_estimate;
    std::cout << "analytic_estimate=" << estimate
              << " pi=" << count_at_estimate
              << " delta=" << (count_at_estimate >= n ? "+" : "-") << initial_delta
              << " setup_ms=" << setup_ms
              << " count_core_ms=" << core_ms
              << " count_ms=" << (setup_ms + core_ms) << '\n';

    const auto start = std::chrono::steady_clock::now();
    const u64 result = nth_prime_counting(n);
    const auto stop = std::chrono::steady_clock::now();
    const auto elapsed =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(stop - start);
    std::cout << "result=" << result << " elapsed_ms=" << elapsed.count() << '\n';
}

u64 nth_prime_with_runtime_state(u64 n) {
    if (n == 0) {
        return 0;
    }
    return nth_prime_counting(n);
}

u64 nth_prime(u64 n) {
    reset_runtime_state();
    const u64 result = nth_prime_with_runtime_state(n);
    reset_runtime_state();
    return result;
}

u64 nth_prime_with_count_config(u64 n, CountConfig config) {
    reset_runtime_state();
    ForcedCountConfigScope scope(config);
    const u64 result = nth_prime_with_runtime_state(n);
    reset_runtime_state();
    return result;
}

void profile_nth_prime(u64 n, bool detail = false) {
    reset_runtime_state();
    g_profile.reset(true, detail);
    g_profile_thread.reset();

    const auto start = std::chrono::steady_clock::now();
    const u64 result = nth_prime_with_runtime_state(n);
    const auto stop = std::chrono::steady_clock::now();

    aggregate_profile_thread();
    g_profile.result = result;
    g_profile.total_ms = std::chrono::duration<double, std::milli>(stop - start).count();

    const auto pct = [](double part, double total) {
        return total > 0.0 ? (100.0 * part / total) : 0.0;
    };
    const double phase_ms[] = {
        g_profile.estimate_ms,
        g_profile.count_setup_ms,
        g_profile.count_core_ms,
        g_profile.correction_ms,
    };
    const char* phase_names[] = {
        "analytic_estimate",
        "count_setup",
        "count_core",
        "correction_sieve",
    };
    std::size_t primary_phase = 0;
    for (std::size_t i = 1; i < sizeof(phase_ms) / sizeof(phase_ms[0]); ++i) {
        if (phase_ms[i] > phase_ms[primary_phase]) {
            primary_phase = i;
        }
    }

    std::cout << result << '\n'
              << "mode: profile-isolated\n"
              << "profile_detail: " << (g_profile.detail ? "on" : "off") << '\n'
              << "n: " << n << '\n'
              << "estimate: " << g_profile.estimate << '\n'
              << "pi_estimate: " << g_profile.count_at_estimate << '\n'
              << "rank_delta: " << g_profile.delta << '\n'
              << "config_sieve_limit: " << g_profile.config_sieve_limit << '\n'
              << "config_rank_limit: " << g_profile.config_rank_limit << '\n'
              << "config_phi_x: " << g_profile.config_phi_x << '\n'
              << "config_phi_a: " << g_profile.config_phi_a << '\n'
              << "config_threads: " << g_profile.config_threads << '\n'
              << "config_a_boost: " << g_profile.config_a_boost << '\n'
              << "config_phi_split_slack: " << g_profile.config_phi_split_slack << '\n'
              << "direct_wheel_a: " << NPRIME_DIRECT_WHEEL_A << '\n'
              << "profile_primary_phase: " << phase_names[primary_phase] << '\n'
              << "profile_primary_phase_ms: " << std::fixed << std::setprecision(6)
              << phase_ms[primary_phase] << '\n'
              << "phase_estimate_pct: " << pct(g_profile.estimate_ms, g_profile.total_ms) << '\n'
              << "phase_count_setup_pct: " << pct(g_profile.count_setup_ms, g_profile.total_ms) << '\n'
              << "phase_count_core_pct: " << pct(g_profile.count_core_ms, g_profile.total_ms) << '\n'
              << "phase_correction_pct: " << pct(g_profile.correction_ms, g_profile.total_ms) << '\n'
              << "estimate_ms: " << std::fixed << std::setprecision(6) << g_profile.estimate_ms << '\n'
              << "counter_sieve_ms: " << g_profile.counter_sieve_ms << '\n'
              << "counter_sieve_alloc_ms: " << g_profile.counter_sieve_alloc_ms << '\n'
              << "counter_sieve_mark_ms: " << g_profile.counter_sieve_mark_ms << '\n'
              << "counter_sieve_pi_ms: " << g_profile.counter_sieve_pi_ms << '\n'
              << "counter_sieve_rank_ms: " << g_profile.counter_sieve_rank_ms << '\n'
              << "counter_wheel_ms: " << g_profile.counter_wheel_ms << '\n'
              << "counter_phi_ms: " << g_profile.counter_phi_ms << '\n'
              << "count_setup_ms: " << g_profile.count_setup_ms << '\n'
              << "count_core_ms: " << g_profile.count_core_ms << '\n'
              << "count_prepare_ms: " << g_profile.count_prepare_ms << '\n'
              << "count_base_ms: " << g_profile.count_base_ms << '\n'
              << "count_parallel_ms: " << g_profile.count_parallel_ms << '\n'
              << "count_thread_launch_ms: " << g_profile.count_thread_launch_ms << '\n'
              << "count_join_ms: " << g_profile.count_join_ms << '\n'
              << "count_worker_min_ms: " << g_profile.count_worker_min_ms << '\n'
              << "count_worker_avg_ms: " << g_profile.count_worker_avg_ms << '\n'
              << "count_worker_max_ms: " << g_profile.count_worker_max_ms << '\n'
              << "count_worker_top_sum_ms: " << g_profile.count_worker_top_sum_ms << '\n'
              << "count_worker_top_avg_ms: " << g_profile.count_worker_top_avg_ms << '\n'
              << "count_worker_top_max_ms: " << g_profile.count_worker_top_max_ms << '\n'
              << "count_worker_phi_split_sum_ms: " << g_profile.count_worker_phi_split_sum_ms << '\n'
              << "count_worker_phi_split_avg_ms: " << g_profile.count_worker_phi_split_avg_ms << '\n'
              << "count_worker_phi_split_max_ms: " << g_profile.count_worker_phi_split_max_ms << '\n'
              << "count_main_phi_split_ms: " << g_profile.count_main_phi_split_ms << '\n'
              << "count_main_top_ms: " << g_profile.count_main_top_ms << '\n'
              << "count_partial_reduce_ms: " << g_profile.count_partial_reduce_ms << '\n'
              << "correction_ms: " << g_profile.correction_ms << '\n'
              << "correction_plan_ms: " << g_profile.correction_plan_ms << '\n'
              << "correction_mark_ms: " << g_profile.correction_mark_ms << '\n'
              << "correction_scan_ms: " << g_profile.correction_scan_ms << '\n'
              << "total_ms: " << g_profile.total_ms << '\n'
              << "lehmer_top_x: " << g_profile.top_x << '\n'
              << "lehmer_a: " << g_profile.top_a << '\n'
              << "lehmer_b: " << g_profile.top_b << '\n'
              << "lehmer_c: " << g_profile.top_c << '\n'
              << "lehmer_top_terms: " << g_profile.top_total << '\n'
              << "lehmer_threads: " << g_profile.top_threads << '\n'
              << "lehmer_chunk: " << g_profile.top_chunk << '\n'
              << "top_iterations: " << g_profile.top_iterations.load(std::memory_order_relaxed) << '\n'
              << "inner_iterations: " << g_profile.inner_iterations.load(std::memory_order_relaxed) << '\n'
              << "top_phi_split_terms: "
              << g_profile.top_phi_split_terms.load(std::memory_order_relaxed) << '\n'
              << "count_calls: " << g_profile.count_calls.load(std::memory_order_relaxed) << '\n'
              << "count_memo: disabled\n"
              << "count_value_direct: " << g_profile.count_value_direct.load(std::memory_order_relaxed) << '\n'
              << "count_value_rank_direct: "
              << g_profile.count_value_rank_direct.load(std::memory_order_relaxed) << '\n'
              << "count_value_recurse: " << g_profile.count_value_recurse.load(std::memory_order_relaxed) << '\n'
              << "phi_calls: " << g_profile.phi_calls.load(std::memory_order_relaxed) << '\n'
              << "phi_cache_hits: " << g_profile.phi_cache_hits.load(std::memory_order_relaxed) << '\n'
              << "phi_wheel_hits: " << g_profile.phi_wheel_hits.load(std::memory_order_relaxed) << '\n'
              << "phi_small_count_hits: "
              << g_profile.phi_small_count_hits.load(std::memory_order_relaxed) << '\n'
              << "phi_split_calls: " << g_profile.phi_split_calls.load(std::memory_order_relaxed) << '\n'
              << "phi_split_terms: " << g_profile.phi_split_terms.load(std::memory_order_relaxed) << '\n'
              << "phi_recurrence_calls: "
              << g_profile.phi_recurrence_calls.load(std::memory_order_relaxed) << '\n';
    std::cout << "count_calls_by_depth:";
    for (std::size_t i = 0; i < kProfileDepthBuckets; ++i) {
        if (i + 1 == kProfileDepthBuckets) {
            std::cout << " d>=" << i << '=';
        } else {
            std::cout << " d" << i << '=';
        }
        std::cout << g_profile.count_calls_by_depth[i].load(std::memory_order_relaxed);
    }
    std::cout << '\n';
    std::cout << "count_calls_by_size:";
    for (std::size_t i = 0; i < kProfileSizeBuckets; ++i) {
        std::cout << ' ' << profile_size_bucket_label(i) << '='
                  << g_profile.count_calls_by_size[i].load(std::memory_order_relaxed);
    }
    std::cout << '\n';
    std::cout << "count_calls_by_stage:";
    for (std::size_t i = 0; i < kProfileStageBuckets; ++i) {
        std::cout << ' ' << profile_stage_label(i) << '='
                  << g_profile.count_calls_by_stage[i].load(std::memory_order_relaxed);
    }
    std::cout << '\n';
    std::cout << "count_value_direct_by_stage:";
    for (std::size_t i = 0; i < kProfileStageBuckets; ++i) {
        std::cout << ' ' << profile_stage_label(i) << '='
                  << g_profile.count_value_direct_by_stage[i].load(std::memory_order_relaxed);
    }
    std::cout << '\n';
    std::cout << "count_value_rank_direct_by_stage:";
    for (std::size_t i = 0; i < kProfileStageBuckets; ++i) {
        std::cout << ' ' << profile_stage_label(i) << '='
                  << g_profile.count_value_rank_direct_by_stage[i].load(std::memory_order_relaxed);
    }
    std::cout << '\n';
    std::cout << "count_value_recurse_by_stage:";
    for (std::size_t i = 0; i < kProfileStageBuckets; ++i) {
        std::cout << ' ' << profile_stage_label(i) << '='
                  << g_profile.count_value_recurse_by_stage[i].load(std::memory_order_relaxed);
    }
    std::cout << '\n';
    std::cout << "count_value_outer_recurse_sum_ms_by_stage:";
    for (std::size_t i = 0; i < kProfileStageBuckets; ++i) {
        const double ms = static_cast<double>(
            g_profile.count_value_outer_recurse_ns_by_stage[i].load(std::memory_order_relaxed)) /
            1'000'000.0;
        std::cout << ' ' << profile_stage_label(i) << '=' << ms;
    }
    std::cout << '\n';
    std::cout << "phi_calls_by_stage:";
    for (std::size_t i = 0; i < kProfileStageBuckets; ++i) {
        std::cout << ' ' << profile_stage_label(i) << '='
                  << g_profile.phi_calls_by_stage[i].load(std::memory_order_relaxed);
    }
    std::cout << '\n';
    if (g_profile.detail) {
        std::cout << "phi_calls_by_a:";
        for (std::size_t i = 0; i < kProfilePhiABuckets; ++i) {
            std::cout << ' ' << profile_phi_a_bucket_label(i) << '='
                      << g_profile.phi_calls_by_a[i].load(std::memory_order_relaxed);
        }
        std::cout << '\n';
        std::cout << "phi_calls_by_size:";
        for (std::size_t i = 0; i < kProfileSizeBuckets; ++i) {
            std::cout << ' ' << profile_size_bucket_label(i) << '='
                      << g_profile.phi_calls_by_size[i].load(std::memory_order_relaxed);
        }
        std::cout << '\n';
    }
    std::cout << "correction_direction: " << g_profile.correction_direction << '\n'
              << "correction_chunks: " << g_profile.correction_chunks.load(std::memory_order_relaxed) << '\n'
              << "correction_odds: " << g_profile.correction_odds.load(std::memory_order_relaxed) << '\n'
              << "correction_mark_entries: "
              << g_profile.correction_mark_entries.load(std::memory_order_relaxed) << '\n'
              << "correction_mark_ops: "
              << g_profile.correction_mark_ops.load(std::memory_order_relaxed) << '\n';

    reset_runtime_state();
    g_profile.reset(false);
    g_profile_thread.reset();
}

void profile_nth_prime_with_config(u64 n, CountConfig config, bool detail = false) {
    ForcedCountConfigScope scope(config);
    profile_nth_prime(n, detail);
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
        {49'999'999, 982'451'629},
        {50'000'000, 982'451'653},
        {99'999'999, 2'038'074'739},
        {100'000'000, 2'038'074'743},
        {199'999'999, 4'222'234'727},
        {200'000'000, 4'222'234'741},
        {799'999'999, 18'054'236'951ULL},
        {800'000'000, 18'054'236'957ULL},
        {999'999'999, 22'801'763'477ULL},
        {1'000'000'000, 22'801'763'489ULL},
        {10'000'000'000ULL, 252'097'800'623ULL},
        {20'000'000'000ULL, 518'649'879'439ULL},
        {24'999'999'999ULL, 654'124'187'837ULL},
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
              << "  " << argv0 << " --bench-config N sieve phi_x phi_a threads a_boost split_slack [iterations]\n"
              << "  " << argv0 << " --bench-config-rank N sieve rank phi_x phi_a threads a_boost split_slack [iterations]\n"
              << "  " << argv0 << " --trace N\n"
              << "  " << argv0 << " --profile N\n"
              << "  " << argv0 << " --profile-detail N\n"
              << "  " << argv0 << " --profile-config N sieve phi_x phi_a threads a_boost split_slack\n"
              << "  " << argv0 << " --profile-config-rank N sieve rank phi_x phi_a threads a_boost split_slack\n"
              << "  " << argv0 << " --profile-detail-config N sieve phi_x phi_a threads a_boost split_slack\n"
              << "  " << argv0 << " --profile-detail-config-rank N sieve rank phi_x phi_a threads a_boost split_slack\n"
              << "  " << argv0 << " --self-test\n";
}

CountConfig parse_count_config(char** argv, int offset) {
    const u32 sieve = parse_u32_arg(argv[offset], "sieve");
    return {
        sieve,
        sieve,
        parse_u32_arg(argv[offset + 1], "phi_x"),
        parse_u32_arg(argv[offset + 2], "phi_a"),
        parse_u32_arg(argv[offset + 3], "threads"),
        parse_u32_arg(argv[offset + 4], "a_boost"),
        parse_u32_arg(argv[offset + 5], "split_slack"),
    };
}

CountConfig parse_count_config_rank(char** argv, int offset) {
    return {
        parse_u32_arg(argv[offset], "sieve"),
        parse_u32_arg(argv[offset + 1], "rank"),
        parse_u32_arg(argv[offset + 2], "phi_x"),
        parse_u32_arg(argv[offset + 3], "phi_a"),
        parse_u32_arg(argv[offset + 4], "threads"),
        parse_u32_arg(argv[offset + 5], "a_boost"),
        parse_u32_arg(argv[offset + 6], "split_slack"),
    };
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
                      << "mode: isolated\n"
                      << "iterations: " << iterations << '\n'
                      << "total_ms: " << std::fixed << std::setprecision(3) << elapsed.count() << '\n'
                      << "avg_ms: " << std::fixed << std::setprecision(6)
                      << elapsed.count() / static_cast<double>(iterations) << '\n';
            return 0;
        }

        if (argc >= 9 && std::string(argv[1]) == "--bench-config") {
            const u64 n = parse_u64(argv[2]);
            const CountConfig config = parse_count_config(argv, 3);
            const u64 iterations = argc >= 10 ? parse_u64(argv[9]) : 5;
            if (iterations == 0) {
                throw std::invalid_argument("iterations must be positive");
            }

            u64 result = 0;
            const auto start = std::chrono::steady_clock::now();
            for (u64 i = 0; i < iterations; ++i) {
                u64 query = n;
                benchmark_barrier(query);
                result = nth_prime_with_count_config(query, config);
                benchmark_barrier(result);
            }
            const auto stop = std::chrono::steady_clock::now();
            const auto elapsed =
                std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(stop - start);

            std::cout << result << '\n'
                      << "mode: isolated-config\n"
                      << "iterations: " << iterations << '\n'
                      << "sieve: " << config.sieve_limit << '\n'
                      << "rank: " << config.rank_limit << '\n'
                      << "phi_x: " << config.phi_x << '\n'
                      << "phi_a: " << config.phi_a << '\n'
                      << "threads: " << config.threads << '\n'
                      << "a_boost: " << config.a_boost << '\n'
                      << "split_slack: " << config.phi_split_slack << '\n'
                      << "total_ms: " << std::fixed << std::setprecision(3) << elapsed.count() << '\n'
                      << "avg_ms: " << std::fixed << std::setprecision(6)
                      << elapsed.count() / static_cast<double>(iterations) << '\n';
            return 0;
        }

        if (argc >= 10 && std::string(argv[1]) == "--bench-config-rank") {
            const u64 n = parse_u64(argv[2]);
            const CountConfig config = parse_count_config_rank(argv, 3);
            const u64 iterations = argc >= 11 ? parse_u64(argv[10]) : 5;
            if (iterations == 0) {
                throw std::invalid_argument("iterations must be positive");
            }

            u64 result = 0;
            const auto start = std::chrono::steady_clock::now();
            for (u64 i = 0; i < iterations; ++i) {
                u64 query = n;
                benchmark_barrier(query);
                result = nth_prime_with_count_config(query, config);
                benchmark_barrier(result);
            }
            const auto stop = std::chrono::steady_clock::now();
            const auto elapsed =
                std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(stop - start);

            std::cout << result << '\n'
                      << "mode: isolated-config-rank\n"
                      << "iterations: " << iterations << '\n'
                      << "sieve: " << config.sieve_limit << '\n'
                      << "rank: " << config.rank_limit << '\n'
                      << "phi_x: " << config.phi_x << '\n'
                      << "phi_a: " << config.phi_a << '\n'
                      << "threads: " << config.threads << '\n'
                      << "a_boost: " << config.a_boost << '\n'
                      << "split_slack: " << config.phi_split_slack << '\n'
                      << "total_ms: " << std::fixed << std::setprecision(3) << elapsed.count() << '\n'
                      << "avg_ms: " << std::fixed << std::setprecision(6)
                      << elapsed.count() / static_cast<double>(iterations) << '\n';
            return 0;
        }

        if (argc == 3 && std::string(argv[1]) == "--trace") {
            trace_nth_prime(parse_u64(argv[2]));
            return 0;
        }

        if (argc == 3 && std::string(argv[1]) == "--profile") {
            profile_nth_prime(parse_u64(argv[2]));
            return 0;
        }

        if (argc == 3 && std::string(argv[1]) == "--profile-detail") {
            profile_nth_prime(parse_u64(argv[2]), true);
            return 0;
        }

        if (argc == 9 && std::string(argv[1]) == "--profile-config") {
            profile_nth_prime_with_config(parse_u64(argv[2]), parse_count_config(argv, 3));
            return 0;
        }

        if (argc == 10 && std::string(argv[1]) == "--profile-config-rank") {
            profile_nth_prime_with_config(parse_u64(argv[2]), parse_count_config_rank(argv, 3));
            return 0;
        }

        if (argc == 9 && std::string(argv[1]) == "--profile-detail-config") {
            profile_nth_prime_with_config(parse_u64(argv[2]), parse_count_config(argv, 3), true);
            return 0;
        }

        if (argc == 10 && std::string(argv[1]) == "--profile-detail-config-rank") {
            profile_nth_prime_with_config(parse_u64(argv[2]), parse_count_config_rank(argv, 3), true);
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
