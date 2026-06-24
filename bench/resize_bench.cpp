/*
 * [BENCH] Resize-strategy performance comparison.
 *
 * Compares the two resize implementations available in this branch:
 *   - resize_inplace()   : single-pass in-place transcription (production path)
 *   - resize_enumerate() : enumerate() + per-element re-insert (naive analog of
 *                          Rsqf::resize, made counter-aware)
 *
 * For each quotient size q in a sweep:
 *   1. Build a Bqf(q, c, k, z) and fill it with random k-mer hashes (+ random
 *      counts) up to its 95% load threshold.
 *   2. Copy that filled state and time each strategy resizing q -> q+1 on
 *      identical input (R repeats, report the minimum = least-noisy estimate).
 *   3. Assert both strategies produce a byte-identical resized filter.
 *
 * Memory awareness (for clusters / large q):
 *   The enumerate strategy materialises a std::map of all elements (~64 B/node),
 *   which dominates peak RAM. Pass a RAM budget (GB) and the sweep stops before
 *   the estimated peak would exceed it. Peak RSS (VmHWM) is reported per row.
 *   Repeats auto-drop to 1 above BIG_N elements (a single run is already long).
 *
 * Output: CSV on stdout ->
 *   q,n_elements,remainder_size,inplace_ms,enumerate_ms,est_peak_gb,peak_rss_gb,repeats_used
 *
 * Usage: ./resize_bench [q_min] [q_max] [repeats] [seed] [ram_budget_gb]
 *   defaults: q_min=14 q_max=31 repeats=3 seed=42 ram_budget_gb=0 (unlimited)
 *   fixed:    c=5  k=32  z=12   (=> s=20, 40-bit hash space, counts in [1,31])
 *
 * On 256 GB the largest feasible point is q=31 (~2.0e9 elements, ~157 GB peak).
 */

#include "abstract_bqf.hpp"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using clk = std::chrono::steady_clock;

// Fixed filter parameters (kept identical to the local run for comparability).
static constexpr uint64_t C_BITS = 5;    // counter bits  -> counts clamp at 31
static constexpr uint64_t K      = 32;   // k-mer size
static constexpr uint64_t Z      = 12;   // s = k - z = 20  -> 2s = 40-bit hashes
static constexpr uint64_t HASH_BITS = 2 * (K - Z);         // = 40

// Per-element cost of the std::map built by resize_enumerate (libstdc++ node:
// 32 B base + 16 B pair + ~16 B malloc overhead). Used only for the estimate.
static constexpr double MAP_NODE_BYTES = 64.0;
// Above this element count, a single resize already takes ~minutes; one repeat
// is enough for a stable estimate, so we stop wasting hours repeating.
static constexpr double BIG_N = 5e8;

static double ms_since(clk::time_point t0) {
    return std::chrono::duration<double, std::milli>(clk::now() - t0).count();
}

// Peak resident set size (high-water mark) so far, in GB. Monotonic over the
// process lifetime, so the value after sweeping q reflects that q's peak.
static double peak_rss_gb() {
    std::ifstream st("/proc/self/status");
    std::string key;
    while (st >> key) {
        if (key == "VmHWM:") {
            double kb;
            st >> kb;
            return kb / (1024.0 * 1024.0);
        }
        st.ignore(1 << 20, '\n');
    }
    return 0.0;
}

// Analytic peak-memory estimate (GB) for the enumerate strategy at quotient q.
// Sums: the element map + old filter + new filter + retained base + retained
// in-place result. The map term dominates for large q.
static double est_peak_gb(uint64_t q) {
    const double slots   = static_cast<double>(1ULL << q);
    const double n       = 0.95 * slots;
    const double wpb_q   = (HASH_BITS + C_BITS - q) + MET_UNIT;   // words/block at q
    const double wpb_q1  = (HASH_BITS + C_BITS - (q + 1)) + MET_UNIT;
    const double filt_q  = slots       * wpb_q  / BLOCK_SIZE * sizeof(uint64_t);
    const double filt_q1 = slots * 2.0 * wpb_q1 / BLOCK_SIZE * sizeof(uint64_t);
    const double map_b   = n * MAP_NODE_BYTES;
    // map + filter(b) + new_filter + base + last_inplace(q+1)
    const double bytes = map_b + filt_q + filt_q1 + filt_q + filt_q1;
    return bytes / 1e9;
}

// Two Bqf states are equal iff their geometry and packed filter match.
static bool same_filter(const Bqf& a, const Bqf& b) {
    return a.quotient_size  == b.quotient_size  &&
           a.remainder_size == b.remainder_size &&
           a.elements_inside == b.elements_inside &&
           a.filter == b.filter;
}

int main(int argc, char** argv) {
    const uint64_t q_min     = (argc > 1) ? std::stoul(argv[1]) : 14;
    const uint64_t q_max     = (argc > 2) ? std::stoul(argv[2]) : 31;
    const int      repeats   = (argc > 3) ? std::stoi(argv[3])  : 3;
    const uint64_t seed      = (argc > 4) ? std::stoul(argv[4]) : 42;
    const double   ram_gb    = (argc > 5) ? std::stod(argv[5])  : 0.0;   // 0 = unlimited

    const uint64_t hash_mask =
        (HASH_BITS >= 64) ? ~0ULL : ((1ULL << HASH_BITS) - 1);

    std::cerr << "# resize bench: q=[" << q_min << ".." << q_max << "] "
              << "c=" << C_BITS << " k=" << K << " z=" << Z
              << " repeats=" << repeats << " seed=" << seed
              << " ram_budget_gb=" << ram_gb << "\n";
    std::cout << "q,n_elements,remainder_size,inplace_ms,enumerate_ms,"
                 "est_peak_gb,peak_rss_gb,repeats_used\n";

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> hash_dist(0, hash_mask);
    std::uniform_int_distribution<uint64_t> count_dist(1, (1ULL << C_BITS) - 1);

    for (uint64_t q = q_min; q <= q_max; ++q) {
        const double est = est_peak_gb(q);
        if (ram_gb > 0.0 && est > ram_gb) {
            std::cerr << "# stopping before q=" << q << ": estimated peak "
                      << est << " GB exceeds budget " << ram_gb << " GB\n";
            break;
        }

        // ---- build a filter filled to just under the resize threshold ----
        Bqf base(q, C_BITS, K, Z, CountMode::ExactCount, false);
        const uint64_t target = base.size_limit - 1;   // next insert would resize

        auto t_fill = clk::now();
        while (base.elements_inside < target) {
            const uint64_t h = hash_dist(rng) & hash_mask;
            base.insert(h, count_dist(rng));           // dups just bump a counter
        }
        const uint64_t n_elements     = base.elements_inside;
        const uint64_t remainder_size = base.remainder_size;
        const int rused = (static_cast<double>(n_elements) > BIG_N)
                              ? std::min(repeats, 1) : repeats;
        std::cerr << "  q=" << q << " filled n=" << n_elements
                  << " in " << ms_since(t_fill) / 1000.0 << "s"
                  << " (est_peak=" << est << "GB, repeats=" << rused << ")\n";

        // ---- time both strategies on identical copies, take the min ----
        double best_inplace = 1e300, best_enum = 1e300;
        Bqf last_inplace, last_enum;

        for (int r = 0; r < rused; ++r) {
            Bqf a = base;                              // fresh copy at size q
            auto t0 = clk::now();
            a.resize_inplace(1);
            best_inplace = std::min(best_inplace, ms_since(t0));
            last_inplace = std::move(a);

            Bqf b = base;
            t0 = clk::now();
            b.resize_enumerate(1);
            best_enum = std::min(best_enum, ms_since(t0));
            last_enum = std::move(b);
        }

        // ---- correctness: both strategies must agree exactly ----
        if (!same_filter(last_inplace, last_enum)) {
            std::cerr << "MISMATCH at q=" << q
                      << " : resize strategies produced different filters!\n";
            return EXIT_FAILURE;
        }

        std::cout << q << "," << n_elements << "," << remainder_size << ","
                  << best_inplace << "," << best_enum << ","
                  << est << "," << peak_rss_gb() << "," << rused
                  << "\n" << std::flush;                // flush: partial-run safe
        std::cerr << "  q=" << q << " n=" << n_elements
                  << " inplace=" << best_inplace << "ms"
                  << " enumerate=" << best_enum << "ms"
                  << "  (x" << (best_enum / best_inplace) << ")"
                  << "  peak_rss=" << peak_rss_gb() << "GB\n";
    }

    return EXIT_SUCCESS;
}
