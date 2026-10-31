import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from xtal import (CU_KA1, bragg_two_theta, cubic_multiplicity, d_spacing, d_spacing_cubic, f0_cromer_mann,
                  is_allowed, powder_peaks, read_cif, structure_factor)
from xtal.extinction import analytic_F_over_f

CIF = {k: os.path.join(ROOT, "data", v) for k, v in
       {"Si": "Si_COD_9011998.cif", "Al": "Al_COD_9012002.cif", "Fe": "Fe_COD_9008536.cif"}.items()}
KIND = {"Si": "diamond", "Al": "fcc", "Fe": "bcc"}
A = {"Si": 5.430941, "Al": 4.04961, "Fe": 2.8665}
NATOMS = {"Si": 8, "Al": 4, "Fe": 2}


@pytest.fixture(scope="module")
def xtals():
    return {k: read_cif(v, k) for k, v in CIF.items()}


def F(x, hkl):
    d = d_spacing(hkl, x.cell)
    return structure_factor(x, [hkl], [1 / (2 * d)])[0]


# ---------- CIF parsing ----------
@pytest.mark.parametrize("el", ["Si", "Al", "Fe"])
def test_cif_cell_and_atom_count(xtals, el):
    x = xtals[el]
    assert x.a == pytest.approx(A[el])
    assert x.n_atoms == NATOMS[el]


@pytest.mark.parametrize("el", ["Si", "Al", "Fe"])
def test_cif_matches_ase(xtals, el):
    from ase.io import read
    at = read(CIF[el])
    assert len(at) == xtals[el].n_atoms
    assert np.allclose(at.cell.cellpar(), xtals[el].cell)


# ---------- geometry ----------
def test_d_spacing_cubic_values():
    assert d_spacing_cubic((1, 1, 1), 5.430941) == pytest.approx(5.430941 / np.sqrt(3))
    assert d_spacing_cubic((1, 1, 0), 2.8665) == pytest.approx(2.026922, abs=1e-6)  # 2.8665/sqrt(2)
    assert d_spacing_cubic((2, 0, 0), 4.04961) == pytest.approx(2.024805, abs=1e-6)


@pytest.mark.parametrize("hkl", [(1, 0, 0), (1, 1, 1), (3, 1, 1), (4, 2, 2), (5, 3, 1)])
def test_general_metric_equals_cubic(hkl):
    cell = (5.430941, 5.430941, 5.430941, 90, 90, 90)
    assert d_spacing(hkl, cell) == pytest.approx(d_spacing_cubic(hkl, 5.430941), rel=1e-12)


def test_general_metric_hexagonal():
    # hexagonal: 1/d^2 = 4/3 (h^2+hk+k^2)/a^2 + l^2/c^2
    a, c = 3.0, 5.0
    h, k, l = 1, 1, 2
    expect = 1 / np.sqrt(4 / 3 * (h * h + h * k + k * k) / a ** 2 + l * l / c ** 2)
    assert d_spacing((h, k, l), (a, a, c, 90, 90, 120)) == pytest.approx(expect, rel=1e-12)


def test_bragg_angle_si_111():
    d = 5.430941 / np.sqrt(3)
    tt = bragg_two_theta(d, CU_KA1)
    assert np.sin(np.radians(tt / 2)) * 2 * d == pytest.approx(CU_KA1, rel=1e-12)
    assert tt == pytest.approx(28.4424, abs=1e-3)


def test_bragg_inaccessible_is_nan():
    assert np.isnan(bragg_two_theta(0.5, CU_KA1))


# ---------- form factor ----------
@pytest.mark.parametrize("el,Z", [("Al", 13), ("Si", 14), ("Fe", 26)])
def test_form_factor_at_zero_equals_Z(el, Z):
    assert f0_cromer_mann(el, 0.0) == pytest.approx(Z, abs=0.01)


def test_form_factor_decreases():
    s = np.linspace(0, 1.5, 50)
    assert np.all(np.diff(f0_cromer_mann("Fe", s)) < 0)


# ---------- extinction rules ----------
@pytest.mark.parametrize("hkl", [(1, 0, 0), (1, 1, 0), (2, 1, 0), (2, 1, 1)])
def test_al_fcc_forbidden(xtals, hkl):
    assert not is_allowed(hkl, "fcc")
    assert abs(F(xtals["Al"], hkl)) < 1e-10


@pytest.mark.parametrize("hkl", [(1, 0, 0), (1, 1, 1), (2, 1, 0), (2, 2, 1)])
def test_fe_bcc_forbidden(xtals, hkl):
    assert not is_allowed(hkl, "bcc")
    assert abs(F(xtals["Fe"], hkl)) < 1e-10


@pytest.mark.parametrize("hkl", [(2, 0, 0), (2, 2, 2), (4, 2, 0), (1, 0, 0), (1, 1, 0)])
def test_si_diamond_forbidden(xtals, hkl):
    assert not is_allowed(hkl, "diamond")
    assert abs(F(xtals["Si"], hkl)) < 1e-10


@pytest.mark.parametrize("el", ["Si", "Al", "Fe"])
def test_numeric_F_matches_analytic_for_all_hkl(xtals, el):
    x = xtals[el]
    for h in range(0, 6):
        for k in range(0, 6):
            for l in range(0, 6):
                if (h, k, l) == (0, 0, 0):
                    continue
                d = d_spacing((h, k, l), x.cell)
                s = 1 / (2 * d)
                f = f0_cromer_mann(el, s)
                Fnum = structure_factor(x, [(h, k, l)], [s])[0]
                # magnitude test (phase depends on origin choice)
                assert abs(Fnum) == pytest.approx(abs(analytic_F_over_f((h, k, l), KIND[el])) * f, abs=1e-8)
                assert (abs(Fnum) > 1e-8) == is_allowed((h, k, l), KIND[el])


def test_si_magnitudes():
    # |F|/f: all-odd 4√2, h+k+l=4n even 8
    assert abs(analytic_F_over_f((1, 1, 1), "diamond")) == pytest.approx(4 * np.sqrt(2))
    assert abs(analytic_F_over_f((2, 2, 0), "diamond")) == pytest.approx(8)


# ---------- multiplicity ----------
@pytest.mark.parametrize("hkl,m", [((1, 0, 0), 6), ((1, 1, 0), 12), ((1, 1, 1), 8), ((2, 1, 0), 24),
                                   ((2, 1, 1), 24), ((3, 2, 1), 48)])
def test_cubic_multiplicity(hkl, m):
    assert cubic_multiplicity(hkl) == m


def test_enumerated_multiplicity_matches_cubic(xtals):
    for el, x in xtals.items():
        for p in powder_peaks(x, CU_KA1, (10, 120)):
            assert p.multiplicity == sum(cubic_multiplicity(f) for f in p.families)


# ---------- cross-validation vs pymatgen ----------
@pytest.mark.parametrize("el", ["Si", "Al", "Fe"])
def test_two_theta_agrees_with_pymatgen(xtals, el):
    from xtal.compare import compare
    rows, own, pat = compare(xtals[el], CIF[el], CU_KA1, (10, 120))
    assert len(own) == len(pat.x)
    for r in rows:
        assert abs(r["dtt"]) < 1e-6


@pytest.mark.parametrize("el", ["Si", "Al", "Fe"])
def test_intensity_close_to_pymatgen(xtals, el):
    """Different form-factor tables -> small differences; require < 2 (on 0–100 scale)."""
    from xtal.compare import compare
    rows, _, _ = compare(xtals[el], CIF[el], CU_KA1, (10, 120))
    assert max(abs(r["dI"]) for r in rows) < 2.0


@pytest.mark.parametrize("el", ["Si", "Al", "Fe"])
def test_intensity_identical_with_same_form_factor(xtals, el):
    """With pymatgen's own f(s) plugged in, the algorithms must agree to round-off."""
    from xtal.compare import compare, pymatgen_form_factor
    rows, _, _ = compare(xtals[el], CIF[el], CU_KA1, (10, 120), form_factor=pymatgen_form_factor)
    assert max(abs(r["dI"]) for r in rows) < 1e-9


def test_expected_peak_lists(xtals):
    first = {el: [p.families[0] for p in powder_peaks(x, CU_KA1, (10, 120))][:4] for el, x in xtals.items()}
    assert first["Si"] == [(1, 1, 1), (2, 2, 0), (3, 1, 1), (4, 0, 0)]
    assert first["Al"] == [(1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1)]
    assert first["Fe"] == [(1, 1, 0), (2, 0, 0), (2, 1, 1), (2, 2, 0)]
