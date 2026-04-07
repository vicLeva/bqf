# Changelog

## [1.1.0] - 2026-04-07

### Breaking changes
- `Bqf_ec` and `Bqf_oom` are removed. Replaced by a single `Bqf` class with a `CountMode` enum (`EC` / `OOM`). Update construction calls accordingly.

### Performance
- Branchless nucleotide encoding via 4-entry LUT (replaces switch)
- Branchless encode/decode using bit ops (replaces switches)
- Incremental reverse complement in `query_sequence`, replacing per-s-mer `revcomp64`
- Precomputed `words_per_block` / `bits_per_block` at construction
- Loop-invariant `smer_mask` and `zp1` hoisted out of `query_sequence` inner loop
- `reserve(k)` call in `decode` / `quick_decoding` paths
- Bit ops in `get_block_id` / `get_shift_in_block`

### Bug fixes
- Fixed undefined behaviour in `bitselectasm` for 0-input (`ctzll(0)`)
- Guard against `quotient_size >= 64` UB in `Rsqf` constructor
- Fixed `result_query` field types (`int`/`float` → `uint64_t`/`double`)
- Fixed `int64_t`/`uint64_t` mismatch in `bqf_oom` `load_from_disk`
- `insert(string)` now throws on s-mer size mismatch instead of silently returning
- Added `file.good()` checks in `load_from_disk` serialization

### Code quality
- API cleanup: const correctness throughout, removed gtest from production code
- Overhauled unit tests: uncommented internal tests, added serialization round-trips, removed dead setup
- `bitrankasm` / `bitselectasm` inlined in header

## [1.0.0] - 2024-02-09

Initial release.
