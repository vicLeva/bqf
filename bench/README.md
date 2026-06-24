# Resize-strategy benchmark

Compares the two BQF resize implementations on this branch:

- **improved resize** (`Bqf::resize_inplace`) — the production path: one pass
  over the old filter, writing remainders + counters straight into the grown
  filter. Linear in element count, needs one extra filter of memory.
- **naive resize** (`Bqf::resize_enumerate`) — counter-aware analog of
  `Rsqf::resize`: dump every `(hash, count)` into a `std::map`, then re-`insert`
  each element. Super-linear in time and needs ~64 B/element of auxiliary memory
  for the map (the dominant cost at scale).

Both produce a byte-identical resized filter; the harness asserts this every run.

## Local run

```bash
cmake -S . -B build_bench -DCMAKE_BUILD_TYPE=Release
cmake --build build_bench --target resize_bench -j$(nproc)

./build_bench/bin/resize_bench 14 22 3 42 > results.csv
python3 bench/plot_resize.py results.csv resize   # -> resize_time.png, resize_memory.png
```

The plot script writes two figures (500 dpi): `<prefix>_time.png` (resize time +
speedup) and `<prefix>_memory.png` (memory required by each variant).

`resize_bench [q_min] [q_max] [repeats] [seed] [ram_budget_gb]`
- Fills a `Bqf(q, c=5, k=32, z=12)` with random k-mer hashes (+random counts) to
  95% load, then times each strategy resizing `q -> q+1`.
- `ram_budget_gb > 0` stops the sweep before the estimated peak would exceed it.
- Repeats auto-drop to 1 above 5e8 elements.

CSV columns:
`q,n_elements,remainder_size,inplace_ms,enumerate_ms,improved_mem_gb,naive_mem_gb,peak_rss_gb,repeats_used`
(`improved_mem_gb`/`naive_mem_gb` are each variant's modelled memory requirement;
`peak_rss_gb` is the measured process high-water mark, for validation.)

## Cluster run (memory-sized)

```bash
# from the repo root, on the cluster node:
bench/run_cluster.sh 256          # RAM_GB [REPEATS] [Q_MIN] [OUTDIR]
```

The bench self-limits `q_max` to fit `RAM_GB - 48` GB of headroom. On 256 GB the
largest point reached is **q = 31** (~2.0e9 elements, ~157 GB peak). The sweep
can take several hours, dominated by the largest `q` (fill + enumerate); run it
under `nohup` or a scheduler. Output (CSV + log + the two PNGs) lands in `OUTDIR`
(default `bench_out/`). Rows are flushed as they complete, so a killed run still
leaves a usable partial CSV.
