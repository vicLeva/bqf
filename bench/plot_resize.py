#!/usr/bin/env python3
"""[BENCH] Plot resize-strategy comparison from resize_bench CSV.

Usage:
    ./bin/resize_bench > results.csv
    python3 bench/plot_resize.py results.csv [out.png]

CSV columns (extra columns are tolerated):
    q,n_elements,remainder_size,inplace_ms,enumerate_ms,
    [est_peak_gb,peak_rss_gb,repeats_used]
"""
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: plot_resize.py <results.csv> [out.png]")
    csv_path = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else "resize_compare.png"

    n, inplace, enumerate_, peak = [], [], [], []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            n.append(int(row["n_elements"]))
            inplace.append(float(row["inplace_ms"]))
            enumerate_.append(float(row["enumerate_ms"]))
            # peak RSS is optional (cluster runs only); fall back to estimate
            val = row.get("peak_rss_gb") or row.get("est_peak_gb")
            peak.append(float(val) if val else None)

    have_mem = all(p is not None for p in peak)
    ncols = 3 if have_mem else 2
    fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 5))
    ax = axes[0]
    ax2 = axes[1]

    # Panel 1: absolute time, log-log
    ax.plot(n, enumerate_, "o-", color="#d1495b",
            label="enumerate + reinsert\n(naive, counter-aware)")
    ax.plot(n, inplace, "s-", color="#2e86ab",
            label="in-place transcription\n(production)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("elements in filter before resize")
    ax.set_ylabel("resize time (ms)")
    ax.set_title("BQF resize: time vs. filter occupancy")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend()

    # Panel 2: speedup factor
    speedup = [e / i for e, i in zip(enumerate_, inplace)]
    ax2.plot(n, speedup, "d-", color="#3c887e")
    ax2.set_xscale("log")
    ax2.set_xlabel("elements in filter before resize")
    ax2.set_ylabel("speedup  (enumerate / in-place)")
    ax2.set_title("In-place speedup over naive enumerate")
    ax2.grid(True, which="both", ls=":", alpha=0.5)
    ax2.axhline(1.0, color="gray", ls="--", lw=0.8)

    # Panel 3 (cluster runs): peak resident memory
    if have_mem:
        ax3 = axes[2]
        ax3.plot(n, peak, "^-", color="#9b59b6")
        ax3.set_xscale("log")
        ax3.set_xlabel("elements in filter before resize")
        ax3.set_ylabel("peak RSS (GB)")
        ax3.set_title("Process peak memory\n(dominated by enumerate map)")
        ax3.grid(True, which="both", ls=":", alpha=0.5)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
