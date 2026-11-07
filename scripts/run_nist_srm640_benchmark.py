"""Increment 2: benchmark COD Si / project Bragg calc against NIST SRM 640e & 640f.

Teaching simulation — we do NOT measure SRM powder. We:
  (a) compare COD Si lattice parameter to NIST certified a;
  (b) recompute Table 1 2θ from NIST certified a + certificate λ and compare to
      NIST information values (should match to ~0.001°);
  (c) illustrate sensitivity by computing 2θ with COD a + project λ=1.5406 Å
      vs NIST 640e Table 1 (different a and slightly different λ — not a claim
      that COD equals SRM).

Usage:  python scripts/run_nist_srm640_benchmark.py
Outputs:
  results/nist_srm640_benchmark.md
  results/nist_srm640e_peaks.csv
  results/nist_srm640f_peaks.csv
  results/nist_srm640_lattice.csv
  figures/Si_nist_srm640e_peak_positions.png
"""
from __future__ import annotations

import csv
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from xtal import CU_KA1, cubic_bragg_two_theta, read_cif

NIST_DIR = os.path.join(ROOT, "data", "nist")
CIF_SI = os.path.join(ROOT, "data", "Si_COD_9011998.cif")


def load_srm(srm_id: str) -> dict:
    with open(os.path.join(NIST_DIR, f"SRM_{srm_id}.json")) as fh:
        return json.load(fh)


def ppm(delta: float, ref: float) -> float:
    return delta / ref * 1.0e6


def peak_rows(reflections, a, wavelength, label_a, label_lam):
    rows = []
    for r in reflections:
        hkl = (r["h"], r["k"], r["l"])
        nist_tt = float(r["two_theta_deg"])
        calc = float(cubic_bragg_two_theta(hkl, a, wavelength))
        rows.append({
            "h": r["h"], "k": r["k"], "l": r["l"],
            "hkl": f"{r['h']}{r['k']}{r['l']}",
            "nist_two_theta_deg": nist_tt,
            "calc_two_theta_deg": calc,
            "delta_two_theta_deg": calc - nist_tt,
            "a_Angstrom": a,
            "wavelength_Angstrom": wavelength,
            "a_source": label_a,
            "wavelength_source": label_lam,
        })
    return rows


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def make_figure(rows_certified, rows_cod, out_path, a_e, lam_nist, a_cod):
    """Two-panel stick comparison: NIST Table 1 vs recomputed; COD+project λ vs NIST."""
    fig, axes = plt.subplots(2, 1, figsize=(10, 7.2), sharex=False)

    # Panel 1: NIST 640e Table 1 sticks vs recomputed from certified a + certificate λ
    ax = axes[0]
    nist_tt = [r["nist_two_theta_deg"] for r in rows_certified]
    calc_tt = [r["calc_two_theta_deg"] for r in rows_certified]
    ax.vlines(nist_tt, 0, 1.0, color="C3", lw=2.0, label="NIST 640e Table 1 (information values)")
    ax.vlines(calc_tt, 0, -1.0, color="C0", lw=2.0, label="recomputed from certified a + certificate λ")
    ax.axhline(0, color="k", lw=0.6)
    for r in rows_certified:
        ax.annotate(r["hkl"], (r["nist_two_theta_deg"], 1.05), ha="center", fontsize=7, rotation=90)
    ax.set_ylim(-1.35, 1.55)
    ax.set_xlim(20, 145)
    ax.set_ylabel("stick (arb.)")
    ax.set_title(
        f"Si NIST SRM 640e — Table 1 vs Bragg recomputation "
        f"(a = {a_e:.6f} Å, λ = {lam_nist:.7f} Å)\n"
        "TEACHING BENCHMARK — transcription of NIST certificate; NOT our measurement"
    )
    ax.legend(loc="lower right", fontsize=8)
    ax.text(0.01, 0.02, "Information values given to 0.001°; recomputation should agree to ~0.001°",
            transform=ax.transAxes, fontsize=8, color="gray", va="bottom")

    # Panel 2: COD a + project λ vs NIST 640e Table 1 (sensitivity illustration)
    ax = axes[1]
    nist_tt = [r["nist_two_theta_deg"] for r in rows_cod]
    calc_tt = [r["calc_two_theta_deg"] for r in rows_cod]
    ax.vlines(nist_tt, 0, 1.0, color="C3", lw=2.0, label="NIST 640e Table 1")
    ax.vlines(calc_tt, 0, -1.0, color="C2", lw=2.0,
              label=f"COD a={a_cod:.6f} Å + project λ={CU_KA1} Å")
    ax.axhline(0, color="k", lw=0.6)
    for r in rows_cod:
        ax.annotate(r["hkl"], (r["nist_two_theta_deg"], 1.05), ha="center", fontsize=7, rotation=90)
    ax.set_ylim(-1.35, 1.55)
    ax.set_xlim(20, 145)
    ax.set_xlabel("2θ (degrees)")
    ax.set_ylabel("stick (arb.)")
    ax.set_title(
        "Sensitivity illustration: COD Si a + project λ vs NIST 640e Table 1\n"
        "(different a and slightly different λ — NOT a claim that COD equals SRM)"
    )
    ax.legend(loc="lower right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    fig_dir = os.path.join(ROOT, "figures")
    res_dir = os.path.join(ROOT, "results")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(res_dir, exist_ok=True)

    srm_e = load_srm("640e")
    srm_f = load_srm("640f")
    a_e = float(srm_e["certified_lattice_parameter"]["value_Angstrom"])
    a_f = float(srm_f["certified_lattice_parameter"]["value_Angstrom"])
    lam_nist = float(srm_e["table1_wavelength"]["value_Angstrom"])
    assert float(srm_f["table1_wavelength"]["value_Angstrom"]) == lam_nist

    xt = read_cif(CIF_SI, "Si")
    a_cod = float(xt.a)

    # (a) lattice comparison
    lattice_rows = [
        {
            "source": "COD 9011998 (Hom, Kiszenick & Post 1975)",
            "a_Angstrom": a_cod,
            "temperature": "298.15 K",
            "delta_vs_640e_Angstrom": a_cod - a_e,
            "delta_vs_640e_ppm": ppm(a_cod - a_e, a_e),
            "delta_vs_640f_Angstrom": a_cod - a_f,
            "delta_vs_640f_ppm": ppm(a_cod - a_f, a_f),
        },
        {
            "source": "NIST SRM 640e certified (22.5 °C)",
            "a_Angstrom": a_e,
            "temperature": "22.5 °C",
            "delta_vs_640e_Angstrom": 0.0,
            "delta_vs_640e_ppm": 0.0,
            "delta_vs_640f_Angstrom": a_e - a_f,
            "delta_vs_640f_ppm": ppm(a_e - a_f, a_f),
        },
        {
            "source": "NIST SRM 640f certified (22.5 °C)",
            "a_Angstrom": a_f,
            "temperature": "22.5 °C",
            "delta_vs_640e_Angstrom": a_f - a_e,
            "delta_vs_640e_ppm": ppm(a_f - a_e, a_e),
            "delta_vs_640f_Angstrom": 0.0,
            "delta_vs_640f_ppm": 0.0,
        },
    ]
    write_csv(
        os.path.join(res_dir, "nist_srm640_lattice.csv"),
        lattice_rows,
        list(lattice_rows[0].keys()),
    )

    # (b) recompute Table 1 from certified a + certificate λ
    rows_e = peak_rows(
        srm_e["table1_peak_positions"]["reflections"],
        a_e, lam_nist, "NIST SRM 640e certified a", "NIST certificate Table 1 λ",
    )
    rows_f = peak_rows(
        srm_f["table1_peak_positions"]["reflections"],
        a_f, lam_nist, "NIST SRM 640f certified a", "NIST certificate Table 1 λ",
    )
    peak_fields = list(rows_e[0].keys())
    write_csv(os.path.join(res_dir, "nist_srm640e_peaks.csv"), rows_e, peak_fields)
    write_csv(os.path.join(res_dir, "nist_srm640f_peaks.csv"), rows_f, peak_fields)

    # (c) COD a + project λ vs NIST 640e Table 1
    rows_cod = peak_rows(
        srm_e["table1_peak_positions"]["reflections"],
        a_cod, CU_KA1, "COD 9011998 a", "project CU_KA1=1.5406 Å",
    )
    write_csv(os.path.join(res_dir, "nist_srm640e_cod_vs_table1.csv"), rows_cod, peak_fields)

    max_d_e = max(abs(r["delta_two_theta_deg"]) for r in rows_e)
    max_d_f = max(abs(r["delta_two_theta_deg"]) for r in rows_f)
    max_d_cod = max(abs(r["delta_two_theta_deg"]) for r in rows_cod)

    fig_path = os.path.join(fig_dir, "Si_nist_srm640e_peak_positions.png")
    make_figure(rows_e, rows_cod, fig_path, a_e, lam_nist, a_cod)

    # Markdown report
    md = []
    md.append("# NIST SRM 640 Si lattice / peak-position benchmark (Increment 2)")
    md.append("")
    md.append("> **Teaching simulation + literature/COD vs NIST certified standard comparison.**")
    md.append("> We did **not** measure SRM powder. Peak lists are NIST *information values*")
    md.append("> transcribed from the certificates. No novelty claimed.")
    md.append("")
    md.append("Auto-generated by `scripts/run_nist_srm640_benchmark.py`.")
    md.append("")
    md.append("## Sources")
    md.append("")
    md.append("| SRM | Certificate | Issue date | Local PDF | SHA-256 |")
    md.append("|---|---|---|---|---|")
    for s in (srm_e, srm_f):
        md.append(
            f"| {s['srm_id']} | [{s['certificate_url']}]({s['certificate_url']}) "
            f"| {s['certificate_issue_date']} | `{s['local_pdf']}` | `{s['sha256']}` |"
        )
    md.append("")
    md.append(f"Download date: **2026-11-07** (Asia/Shanghai). See also `data/SOURCES.md`.")
    md.append("")
    md.append("## (a) Lattice parameter comparison")
    md.append("")
    md.append("| Source | a (Å) | T | Δa vs 640e (Å) | Δa vs 640e (ppm) | Δa vs 640f (Å) | Δa vs 640f (ppm) |")
    md.append("|---|---:|---|---:|---:|---:|---:|")
    for r in lattice_rows:
        md.append(
            f"| {r['source']} | {r['a_Angstrom']:.6f} | {r['temperature']} "
            f"| {r['delta_vs_640e_Angstrom']:+.6f} | {r['delta_vs_640e_ppm']:+.2f} "
            f"| {r['delta_vs_640f_Angstrom']:+.6f} | {r['delta_vs_640f_ppm']:+.2f} |"
        )
    md.append("")
    md.append(
        f"COD a = **{a_cod:.6f} Å** (298.15 K) is **{abs(ppm(a_cod - a_e, a_e)):.2f} ppm** "
        f"smaller than NIST 640e certified a = **{a_e:.6f} Å** (22.5 °C), and "
        f"**{abs(ppm(a_cod - a_f, a_f)):.2f} ppm** smaller than 640f a = **{a_f:.6f} Å**. "
        "Temperatures differ slightly (298.15 K ≈ 25 °C vs 22.5 °C); no thermal correction applied here."
    )
    md.append("")
    md.append("## (b) Recompute NIST Table 1 from certified a + certificate λ")
    md.append("")
    md.append(
        f"λ (certificate) = **{lam_nist:.7f} Å** (= 0.15405929 nm). "
        "Expectation: |Δ2θ| ≲ 0.001° (Table 1 is quoted to 0.001°)."
    )
    md.append("")
    md.append("### SRM 640e")
    md.append("")
    md.append("| hkl | NIST Table 1 2θ (°) | recomputed 2θ (°) | Δ2θ (°) |")
    md.append("|---|---:|---:|---:|")
    for r in rows_e:
        md.append(
            f"| {r['hkl']} | {r['nist_two_theta_deg']:.3f} "
            f"| {r['calc_two_theta_deg']:.6f} | {r['delta_two_theta_deg']:+.6f} |"
        )
    md.append("")
    md.append(f"**max \\|Δ2θ\\| (640e)** = **{max_d_e:.6f}°**")
    md.append("")
    md.append("### SRM 640f")
    md.append("")
    md.append("| hkl | NIST Table 1 2θ (°) | recomputed 2θ (°) | Δ2θ (°) |")
    md.append("|---|---:|---:|---:|")
    for r in rows_f:
        md.append(
            f"| {r['hkl']} | {r['nist_two_theta_deg']:.3f} "
            f"| {r['calc_two_theta_deg']:.6f} | {r['delta_two_theta_deg']:+.6f} |"
        )
    md.append("")
    md.append(f"**max \\|Δ2θ\\| (640f)** = **{max_d_f:.6f}°**")
    md.append("")
    md.append("## (c) Sensitivity: COD a + project λ vs NIST 640e Table 1")
    md.append("")
    md.append(
        f"Uses COD a = {a_cod:.6f} Å and project λ = {CU_KA1} Å (Cu Kα1 approximate). "
        "This **mixes different a and slightly different λ** relative to the NIST certificate "
        "— teaching illustration of peak-position sensitivity, **not** a claim that COD equals SRM."
    )
    md.append("")
    md.append("| hkl | NIST 640e Table 1 2θ (°) | COD+projectλ 2θ (°) | Δ2θ (°) |")
    md.append("|---|---:|---:|---:|")
    for r in rows_cod:
        md.append(
            f"| {r['hkl']} | {r['nist_two_theta_deg']:.3f} "
            f"| {r['calc_two_theta_deg']:.6f} | {r['delta_two_theta_deg']:+.6f} |"
        )
    md.append("")
    md.append(f"**max \\|Δ2θ\\| (COD+projectλ vs 640e Table 1)** = **{max_d_cod:.6f}°**")
    md.append("")
    md.append("## Outputs")
    md.append("")
    md.append("- `results/nist_srm640_lattice.csv`")
    md.append("- `results/nist_srm640e_peaks.csv`")
    md.append("- `results/nist_srm640f_peaks.csv`")
    md.append("- `results/nist_srm640e_cod_vs_table1.csv`")
    md.append("- `figures/Si_nist_srm640e_peak_positions.png`")
    md.append("")
    md.append("## Limitations")
    md.append("")
    md.append("- No experimental measurement of SRM powder was performed.")
    md.append("- NIST Table 1 entries are *information values* (no stated uncertainty), not certified values.")
    md.append("- Only cubic Bragg geometry (peak positions) is checked here; intensities, Debye–Waller,")
    md.append("  absorption, Kα2, and instrument function are out of scope for Increment 2.")
    md.append("- COD and NIST lattice parameters refer to slightly different temperatures; no CTE correction.")
    md.append("")

    md_path = os.path.join(res_dir, "nist_srm640_benchmark.md")
    with open(md_path, "w") as fh:
        fh.write("\n".join(md) + "\n")

    print(f"Wrote {md_path}")
    print(f"Wrote {fig_path}")
    print(f"COD a={a_cod:.6f} Å; 640e a={a_e:.6f} Å; 640f a={a_f:.6f} Å")
    print(f"Δa COD−640e = {a_cod - a_e:+.6f} Å ({ppm(a_cod - a_e, a_e):+.2f} ppm)")
    print(f"Δa COD−640f = {a_cod - a_f:+.6f} Å ({ppm(a_cod - a_f, a_f):+.2f} ppm)")
    print(f"max|Δ2θ| 640e recomputed = {max_d_e:.6f}°")
    print(f"max|Δ2θ| 640f recomputed = {max_d_f:.6f}°")
    print(f"max|Δ2θ| COD+projectλ vs 640e = {max_d_cod:.6f}°")


if __name__ == "__main__":
    main()
