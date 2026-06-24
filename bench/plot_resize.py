#!/usr/bin/env python3
"""[BENCH] Plot resize-strategy comparison from resize_bench CSV.

Produces two figures:
    <prefix>_time.png   : resize time (log-log) + speedup
    <prefix>_memory.png : memory required by each variant

Usage:
    ./bin/resize_bench > results.csv
    python3 bench/plot_resize.py results.csv [out_prefix]   # default prefix: resize

Memory per variant is computed analytically from q / remainder_size / n_elements,
so this works on any resize_bench CSV (old or new column set).

CSV columns required:
    q,n_elements,remainder_size,inplace_ms,enumerate_ms,...
"""
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DPI = 500
NAIVE = "naive resize"
IMPROVED = "improved resize"
C_NAIVE = "#d1495b"
C_IMPROVED = "#2e86ab"

MAP_NODE_BYTES = 64  # libstdc++ std::map<u64,u64> node (32 base + 16 pair + ~16 malloc)
MET_UNIT = 3         # metadata words per block (offset/occupied/runend)
BLOCK_SIZE = 64      # slots per block


def filt_bytes(q, remainder_size):
    """Bytes of the packed filter at quotient q with the given remainder size."""
    words_per_block = remainder_size + MET_UNIT
    return (1 << q) * words_per_block / BLOCK_SIZE * 8


def improved_mem_gb(q, remainder_size):
    # old filter (q) + new filter (q+1, remainder shrinks by 1)
    return (filt_bytes(q, remainder_size) + filt_bytes(q + 1, remainder_size - 1)) / 1e9


def naive_mem_gb(q, remainder_size, n):
    # improved working set + the n-element std::map
    return improved_mem_gb(q, remainder_size) + n * MAP_NODE_BYTES / 1e9


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: plot_resize.py <results.csv> [out_prefix]")
    csv_path = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else "resize"

    n, naive_t, improved_t = [], [], []
    naive_m, improved_m = [], []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            q = int(row["q"])
            ne = int(row["n_elements"])
            rsize = int(row["remainder_size"])
            n.append(ne)
            improved_t.append(float(row["inplace_ms"]))    # in-place == improved
            naive_t.append(float(row["enumerate_ms"]))     # enumerate == naive
            improved_m.append(improved_mem_gb(q, rsize))
            naive_m.append(naive_mem_gb(q, rsize, ne))

    # ---- figure 1: time + speedup ----
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax.plot(n, naive_t, "o-", color=C_NAIVE, label=NAIVE)
    ax.plot(n, improved_t, "s-", color=C_IMPROVED, label=IMPROVED)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("elements in filter before resize")
    ax.set_ylabel("resize time (ms)")
    ax.set_title("BQF resize: time vs. filter occupancy")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend()

    speedup = [b / a for b, a in zip(naive_t, improved_t)]
    ax2.plot(n, speedup, "d-", color="#3c887e")
    ax2.set_xscale("log")
    ax2.set_xlabel("elements in filter before resize")
    ax2.set_ylabel("speedup  (naive / improved)")
    ax2.set_title("Improved speedup over naive resize")
    ax2.grid(True, which="both", ls=":", alpha=0.5)
    ax2.axhline(1.0, color="gray", ls="--", lw=0.8)

    fig.tight_layout()
    time_path = f"{prefix}_time.png"
    fig.savefig(time_path, dpi=DPI)
    print(f"wrote {time_path}")

    # ---- figure 2: memory required by each variant ----
    figm, axm = plt.subplots(figsize=(7, 5))
    axm.plot(n, naive_m, "o-", color=C_NAIVE, label=NAIVE)
    axm.plot(n, improved_m, "s-", color=C_IMPROVED, label=IMPROVED)
    axm.set_xscale("log")
    axm.set_yscale("log")
    axm.set_xlabel("elements in filter before resize")
    axm.set_ylabel("memory required (GB)")
    axm.set_title("BQF resize: memory required per variant")
    axm.grid(True, which="both", ls=":", alpha=0.5)
    axm.legend()

    figm.tight_layout()
    mem_path = f"{prefix}_memory.png"
    figm.savefig(mem_path, dpi=DPI)
    print(f"wrote {mem_path}")


if __name__ == "__main__":
    main()
