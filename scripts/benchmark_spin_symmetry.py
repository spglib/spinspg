"""Benchmark the public rutile example, with imports outside measured regions."""

# ruff: noqa: E402 -- Set BLAS thread limits before importing NumPy.

import os

for name in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[name] = "1"

import argparse
import cProfile
import json
import platform
import pstats
from importlib.metadata import version
from pathlib import Path
from statistics import median
from time import perf_counter

import numpy as np

import spinspg


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=["spglib", "moyopy"], default="spglib")
    parser.add_argument("--frames", type=int, default=256)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile-dir", type=Path)
    args = parser.parse_args()

    x = 0.695169
    lattice = np.diag([4.87, 4.87, 3.31])
    positions = np.array(
        [
            [0, 0, 0],
            [0.5, 0.5, 0.5],
            [x, x, 0],
            [-x, -x, 0],
            [-x + 0.5, x + 0.5, 0.5],
            [x + 0.5, -x + 0.5, 0.5],
        ]
    )
    numbers = np.array([0, 0, 1, 1, 1, 1])
    fields = np.zeros((args.frames, 6, 3))
    angles = np.linspace(0, 2 * np.pi, args.frames, endpoint=False)
    fields[:, 0, 0] = 2.5 * np.sin(angles)
    fields[:, 0, 2] = 2.5 * np.cos(angles)
    fields[:, 1] = -fields[:, 0]
    cells = [lattice * scale for scale in np.linspace(0.9, 1.1, args.frames)]
    options = dict(symprec=1e-3, backend=args.backend)

    def repeated():
        for moments in fields:
            spinspg.get_spin_symmetry(lattice, positions, numbers, moments, **options)

    def cold():
        for cell, moments in zip(cells, fields):
            spinspg.get_spin_symmetry(cell, positions, numbers, moments, **options)

    start = perf_counter()
    spinspg.get_spin_symmetry(lattice, positions, numbers, fields[0], **options)
    representative = perf_counter() - start
    print(f"Representative search: {representative * 1e3:.3f} ms")
    print(f"Estimated time per {args.frames}-frame pass: {representative * args.frames:.2f} s")
    measurements = {
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            **{name: version(name) for name in ["numpy", "spglib", "moyopy", "spinspg"]},
        },
        "backend": args.backend,
        "frames": args.frames,
        "repeats": args.repeats,
    }
    workloads = {"repeated": repeated, "cold_distinct_cells": cold}
    if hasattr(spinspg, "prepare_spin_symmetry"):
        start = perf_counter()
        prepared = spinspg.prepare_spin_symmetry(lattice, positions, numbers, **options)
        measurements["preparation_seconds"] = perf_counter() - start

        def reuse():
            for moments in fields:
                prepared.get_spin_symmetry(moments)

        workloads["prepared"] = reuse

    for name, run in workloads.items():
        times = []
        for _ in range(args.repeats):
            start = perf_counter()
            run()
            times.append(perf_counter() - start)
        measurements[name] = {"seconds": times, "median_seconds": median(times)}
        print(f"{name}: {median(times):.4f} s ({median(times) / args.frames * 1e3:.3f} ms/frame)")

    # Profiling is a separate pass and never contributes to the wall-clock timings.
    if args.profile_dir is not None:
        args.profile_dir.mkdir(parents=True, exist_ok=True)
        for name, run in workloads.items():
            profile = cProfile.Profile()
            profile.runcall(run)
            profile.dump_stats(str(args.profile_dir / f"{name}.prof"))
            print(f"\n{name} profile:")
            pstats.Stats(profile).strip_dirs().sort_stats("cumulative").print_stats(20)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(measurements, indent=2) + "\n")


if __name__ == "__main__":
    main()
