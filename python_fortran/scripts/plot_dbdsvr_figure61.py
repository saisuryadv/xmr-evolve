#!/usr/bin/env python3
"""Create a LAW 166 Figure 6.1-style DBDSVR/DBDSDC speedup plot."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator


SIZES = (500, 1000, 2000, 3000)
METHODS = ("DBDSDC", "DBDSVR")
EMAIL16_VARIETIES = {
    1: "zero bidiagonal",
    2: "identity",
    3: "arithmetic bidiagonal",
    4: "geometric bidiagonal",
    5: "log-distributed bidiagonal",
    6: "arithmetic tridiagonal, Cholesky-lifted",
    7: "geometric tridiagonal, Cholesky-lifted",
    8: "clustered tridiagonal, Cholesky-lifted",
    9: "arithmetic bidiagonal × √overflow",
    10: "arithmetic bidiagonal × √underflow",
    11: "arithmetic diagonal",
    12: "geometric diagonal",
    13: "clustered diagonal",
    14: "arithmetic diagonal × √overflow",
    15: "arithmetic diagonal × √underflow",
    16: "SPD geometric tridiagonal, Cholesky-lifted",
}

PAPER29_VARIETIES = {
    1: "110 ones spectrum",
    2: "111 uniform, ε apart",
    3: "112 uniform, ε¹ᐟ⁴ apart",
    4: "113 uniform ε → 1",
    5: "114 uniform ε → 1, signed",
    6: "115 geometric ε → 1",
    7: "116 geometric ε → 1, signed",
    8: "117 random spectrum",
    9: "118 clustered at 1",
    10: "119 clustered at ±1",
    11: "120 clustered at ε",
    12: "121 clustered at ±ε",
    13: "200 ABCON0",
    14: "201 ABCON1",
    15: "202 ABCON2",
    16: "203 ABCON3",
    17: "210 random bidiagonal",
    18: "220 GRADP",
    19: "221 GRADM",
    20: "222 Wilkinson+",
    21: "223 Wilkinson−",
    22: "224 Wilkinson W",
    23: "225 double Wilkinson W",
    24: "230 Clement",
    25: "240 GRO0",
    26: "241 GRO1",
    27: "242 GRO2",
    28: "243 GRO3",
    29: "244 bidiagonal Wilkinson+",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--suite", choices=("email16", "paper29"), default="email16"
    )
    parser.add_argument("--compiler-label", default="GNU Fortran 12.2")
    return parser.parse_args()


def load_rows(raw_dir: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(raw_dir.glob("n*_*.csv")):
        with path.open(newline="") as handle:
            for raw in csv.DictReader(handle):
                rows.append(
                    {
                        "method": raw["method"],
                        "type": int(raw["type"]),
                        "n": int(raw["n"]),
                        "range": raw["range"],
                        "ns": int(raw["ns"]),
                        "time_seconds": float(raw["time_seconds"]),
                        "info": int(raw["info"]),
                        "repetitions": int(raw["repetitions"]),
                    }
                )
    return rows


def validate(
    rows: list[dict[str, object]], varieties: dict[int, str]
) -> None:
    expected = len(SIZES) * len(varieties) * len(METHODS)
    if len(rows) != expected:
        raise SystemExit(f"expected {expected} rows, found {len(rows)}")
    keys = {
        (int(row["n"]), int(row["type"]), str(row["method"])) for row in rows
    }
    if len(keys) != expected:
        raise SystemExit("duplicate or missing (N, type, method) rows")
    for row in rows:
        if int(row["n"]) not in SIZES:
            raise SystemExit(f"unexpected N={row['n']}")
        if int(row["type"]) not in varieties:
            raise SystemExit(f"unexpected type={row['type']}")
        if str(row["method"]) not in METHODS:
            raise SystemExit(f"unexpected method={row['method']}")
        if str(row["range"]) != "A":
            raise SystemExit("figure requires RANGE=A only")
        if int(row["info"]) != 0:
            raise SystemExit(f"nonzero INFO in {row}")
        if int(row["ns"]) != int(row["n"]):
            raise SystemExit(f"MFOUND mismatch in {row}")
        if float(row["time_seconds"]) <= 0.0:
            raise SystemExit(f"nonpositive runtime in {row}")


def paired_points(
    rows: list[dict[str, object]], varieties: dict[int, str]
) -> list[dict[str, object]]:
    keyed = {
        (int(row["n"]), int(row["type"]), str(row["method"])): row
        for row in rows
    }
    points: list[dict[str, object]] = []
    point = 0
    for n in SIZES:
        for type_number, variety in varieties.items():
            point += 1
            dc = keyed[n, type_number, "DBDSDC"]
            vr = keyed[n, type_number, "DBDSVR"]
            dc_seconds = float(dc["time_seconds"])
            vr_seconds = float(vr["time_seconds"])
            speedup = dc_seconds / vr_seconds
            points.append(
                {
                    "point": point,
                    "n": n,
                    "type": type_number,
                    "variety": variety,
                    "dbdsdc_seconds": dc_seconds,
                    "dbdsvr_seconds": vr_seconds,
                    "speedup_dbdsdc_over_dbdsvr": speedup,
                    "winner": "DBDSVR" if speedup > 1.0 else "DBDSDC",
                    "dbdsdc_repetitions": int(dc["repetitions"]),
                    "dbdsvr_repetitions": int(vr["repetitions"]),
                }
            )
    return points


def attach_backend_paths(
    points: list[dict[str, object]], audit_dir: Path
) -> None:
    audit_paths = sorted(audit_dir.glob("backend_audit_n*.csv"))
    if not audit_paths:
        return
    keyed: dict[tuple[int, int], int] = {}
    for path in audit_paths:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                keyed[int(row["n"]), int(row["type"])] = int(
                    row["backend_path"]
                )
    if len(keyed) != len(points):
        raise SystemExit(
            f"expected {len(points)} backend-audit rows, found {len(keyed)}"
        )
    for point in points:
        point["dbdsvr_backend_path"] = keyed[
            int(point["n"]), int(point["type"])
        ]


def write_points(points: list[dict[str, object]], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(points[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(points)


def make_summary(
    points: list[dict[str, object]], varieties: dict[int, str], suite: str
) -> dict[str, object]:
    by_size: dict[str, object] = {}
    for n in SIZES:
        selected = [point for point in points if int(point["n"]) == n]
        ratios = [float(point["speedup_dbdsdc_over_dbdsvr"]) for point in selected]
        by_size[str(n)] = {
            "cases": len(selected),
            "dbdsvr_wins": sum(value > 1.0 for value in ratios),
            "median_speedup": statistics.median(ratios),
            "minimum_speedup": min(ratios),
            "maximum_speedup": max(ratios),
            "dbdsdc_fallbacks": sum(
                int(point.get("dbdsvr_backend_path", 1)) != 1
                for point in selected
            ),
        }
    ratios = [float(point["speedup_dbdsdc_over_dbdsvr"]) for point in points]
    return {
        "definition": "speedup = DBDSDC seconds / DBDSVR seconds",
        "suite": suite,
        "range": "A",
        "sizes": list(SIZES),
        "varieties_per_size": len(varieties),
        "points": len(points),
        "info_failures": 0,
        "mfound_mismatches": 0,
        "dbdsvr_wins": sum(value > 1.0 for value in ratios),
        "median_speedup": statistics.median(ratios),
        "minimum_speedup": min(ratios),
        "maximum_speedup": max(ratios),
        "dbdsdc_fallbacks": sum(
            int(point.get("dbdsvr_backend_path", 1)) != 1
            for point in points
        ),
        "by_size": by_size,
    }


def plot(
    points: list[dict[str, object]],
    output: Path,
    varieties: dict[int, str],
    suite: str,
    compiler_label: str,
) -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#292929",
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "xtick.color": "#292929",
            "ytick.color": "#292929",
        }
    )
    fig, ax = plt.subplots(figsize=(11.5, 7.4))
    fig.subplots_adjust(left=0.105, right=0.985, top=0.82, bottom=0.25)
    native = [
        point
        for point in points
        if int(point.get("dbdsvr_backend_path", 1)) == 1
    ]
    fallback = [
        point
        for point in points
        if int(point.get("dbdsvr_backend_path", 1)) != 1
    ]
    ax.scatter(
        [int(point["point"]) for point in native],
        [float(point["speedup_dbdsdc_over_dbdsvr"]) for point in native],
        marker="*",
        s=88,
        color="#1261A0",
        edgecolors="#083B66",
        linewidths=0.35,
        zorder=3,
        label="native DBDSVDMR3 path",
    )
    if fallback:
        ax.scatter(
            [int(point["point"]) for point in fallback],
            [
                float(point["speedup_dbdsdc_over_dbdsvr"])
                for point in fallback
            ],
            marker="X",
            s=76,
            color="#C44E23",
            edgecolors="#713019",
            linewidths=0.45,
            zorder=4,
            label="DBDSVR output-audit fallback to DBDSDC",
        )
    ax.set_yscale("log")
    ax.set_ylim(0.08, 100.0)
    ax.set_xlim(0.25, len(points) + 0.75)
    ax.axhline(1.0, color="#343434", linewidth=1.25, zorder=2)
    cases_per_size = len(varieties)
    for group in range(1, len(SIZES)):
        boundary = group * cases_per_size + 0.5
        ax.axvline(boundary, color="#AAA69D", linewidth=0.9, linestyle=(0, (1, 4)))
    centers = tuple(
        group * cases_per_size + (cases_per_size + 1) / 2
        for group in range(len(SIZES))
    )
    labels = []
    for n in SIZES:
        selected = [point for point in points if int(point["n"]) == n]
        wins = sum(float(point["speedup_dbdsdc_over_dbdsvr"]) > 1.0 for point in selected)
        labels.append(f"N = {n}\nDBDSVR wins {wins}/{cases_per_size}")
    ax.set_xticks(centers, labels)
    ax.tick_params(axis="x", length=0, pad=10)
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=5))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}×"))
    ax.grid(axis="y", which="major", color="#D7D3CB", linewidth=0.8)
    ax.grid(axis="y", which="minor", color="#EEECE7", linewidth=0.45)
    if fallback:
        ax.legend(loc="lower left", frameon=True, fontsize=9.5)
    ax.set_ylabel("Speedup = DBDSDC time / DBDSVR time")
    suite_label = (
        "email matrix varieties 1 → 16"
        if suite == "email16"
        else "29-case DMATGEN reconstruction, IDs 110 → 244"
    )
    ax.set_xlabel(f"Within each size group: {suite_label}", labelpad=12)
    ax.set_title(
        f"DBDSVR speedup over DBDSDC on the {cases_per_size}-case suite\n"
        "Points above 1× favor DBDSVR; points below 1× favor DBDSDC",
        loc="left",
        fontsize=16,
        fontweight="bold",
        pad=16,
    )
    ax.text(
        0.995,
        0.985,
        "RANGE = A · all singular vectors",
        transform=ax.transAxes,
        ha="right",
        va="top",
        color="#5A5752",
        fontsize=10,
    )
    fig.text(
        0.01,
        0.055,
        "c2 · AMD EPYC 9454 · exactly 1 CPU quota · Reference LAPACK/BLAS 3.11.0 "
        f"· {compiler_label} · −O2 · adaptive ≥0.5 s timing per case",
        ha="left",
        va="bottom",
        color="#625F59",
        fontsize=9.2,
    )
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    varieties = (
        EMAIL16_VARIETIES if args.suite == "email16" else PAPER29_VARIETIES
    )
    rows = load_rows(args.raw_dir)
    validate(rows, varieties)
    points = paired_points(rows, varieties)
    attach_backend_paths(points, args.output_dir)
    write_points(points, args.output_dir / "figure61_points.csv")
    summary = make_summary(points, varieties, args.suite)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    figure_name = (
        "dbdsvr_vs_dbdsdc_figure61.png"
        if args.suite == "email16"
        else "dbdsvr_vs_dbdsdc_figure61_paper29.png"
    )
    plot(
        points,
        args.output_dir / figure_name,
        varieties,
        args.suite,
        args.compiler_label,
    )


if __name__ == "__main__":
    main()
