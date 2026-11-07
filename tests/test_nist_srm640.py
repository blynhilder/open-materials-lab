"""Increment 2 tests: NIST SRM 640e/640f certified a and Table 1 peak-position benchmark."""
from __future__ import annotations

import json
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from xtal import CU_KA1, cubic_bragg_two_theta, read_cif

NIST_DIR = os.path.join(ROOT, "data", "nist")
CIF_SI = os.path.join(ROOT, "data", "Si_COD_9011998.cif")

# Exact values from NIST certificates (fetched 2026-11-07); must match JSON transcriptions.
EXPECTED = {
    "640e": {
        "a_A": 5.431179,
        "a_nm": 0.5431179,
        "lam_A": 1.5405929,
        "sha256": "87c21e029cfd1b00ea7a97add759a3d1730b49f629ccd83c15de90adf1d29042",
        "peaks": {
            (1, 1, 1): 28.441,
            (2, 2, 0): 47.300,
            (3, 1, 1): 56.120,
            (4, 0, 0): 69.126,
            (3, 3, 1): 76.372,
            (4, 2, 2): 88.025,
            (5, 1, 1): 94.947,
            (4, 4, 0): 106.701,
            (5, 3, 1): 114.084,
            (6, 2, 0): 127.534,
            (5, 3, 3): 136.880,
        },
    },
    "640f": {
        "a_A": 5.431144,
        "a_nm": 0.5431144,
        "lam_A": 1.5405929,
        "sha256": "219e337b4eaec8b054fe3f568955700e8866411750eab3975c77cd047d8878ff",
        "peaks": {
            (1, 1, 1): 28.441,
            (2, 2, 0): 47.301,
            (3, 1, 1): 56.120,
            (4, 0, 0): 69.127,
            (3, 3, 1): 76.373,
            (4, 2, 2): 88.026,
            (5, 1, 1): 94.947,
            (4, 4, 0): 106.702,
            (5, 3, 1): 114.085,
            (6, 2, 0): 127.535,
            (5, 3, 3): 136.882,
        },
    },
}


def load_srm(srm_id: str) -> dict:
    with open(os.path.join(NIST_DIR, f"SRM_{srm_id}.json")) as fh:
        return json.load(fh)


@pytest.mark.parametrize("srm_id", ["640e", "640f"])
def test_certified_a_matches_json(srm_id):
    data = load_srm(srm_id)
    lp = data["certified_lattice_parameter"]
    exp = EXPECTED[srm_id]
    assert lp["value_Angstrom"] == pytest.approx(exp["a_A"])
    assert lp["value_nm"] == pytest.approx(exp["a_nm"])
    assert lp["uncertainty_nm"] == pytest.approx(0.000008)
    assert lp["uncertainty_Angstrom"] == pytest.approx(0.00008)
    assert data["table1_wavelength"]["value_Angstrom"] == pytest.approx(exp["lam_A"])
    assert data["sha256"] == exp["sha256"]


@pytest.mark.parametrize("srm_id", ["640e", "640f"])
def test_table1_transcription_matches_expected(srm_id):
    data = load_srm(srm_id)
    got = {
        (r["h"], r["k"], r["l"]): r["two_theta_deg"]
        for r in data["table1_peak_positions"]["reflections"]
    }
    assert got == EXPECTED[srm_id]["peaks"]


@pytest.mark.parametrize("srm_id", ["640e", "640f"])
def test_recompute_table1_from_certified_a_and_lambda(srm_id):
    """Bragg recomputation must reproduce NIST Table 1 to within ~0.001°."""
    data = load_srm(srm_id)
    a = data["certified_lattice_parameter"]["value_Angstrom"]
    lam = data["table1_wavelength"]["value_Angstrom"]
    max_abs = 0.0
    for r in data["table1_peak_positions"]["reflections"]:
        hkl = (r["h"], r["k"], r["l"])
        calc = float(cubic_bragg_two_theta(hkl, a, lam))
        delta = calc - float(r["two_theta_deg"])
        max_abs = max(max_abs, abs(delta))
        assert abs(delta) < 0.001, (
            f"{srm_id} {hkl}: calc={calc:.6f} NIST={r['two_theta_deg']:.3f} "
            f"Δ={delta:+.6f}° exceeds 0.001° — diagnose wavelength/a/rounding, do not loosen"
        )
    assert max_abs < 0.001


def test_cod_si_lattice_parameter():
    xt = read_cif(CIF_SI, "Si")
    assert xt.a == pytest.approx(5.430941)


def test_certificate_pdfs_exist_and_sha256():
    import hashlib
    for srm_id in ("640e", "640f"):
        data = load_srm(srm_id)
        pdf = os.path.join(ROOT, data["local_pdf"])
        assert os.path.isfile(pdf), f"missing {pdf}"
        h = hashlib.sha256()
        with open(pdf, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 16), b""):
                h.update(chunk)
        assert h.hexdigest() == data["sha256"]


def test_cubic_bragg_helper_matches_composition():
    from xtal.lattice import bragg_two_theta, d_spacing_cubic
    hkl, a, lam = (1, 1, 1), 5.431179, 1.5405929
    assert cubic_bragg_two_theta(hkl, a, lam) == pytest.approx(
        float(bragg_two_theta(d_spacing_cubic(hkl, a), lam))
    )


def test_project_lambda_constant():
    assert CU_KA1 == 1.5406
