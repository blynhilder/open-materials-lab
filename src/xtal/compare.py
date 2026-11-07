"""与 pymatgen XRDCalculator 交叉验证 (cross-validation). 仅此模块导入 pymatgen."""
from __future__ import annotations

import numpy as np

from .powder import powder_peaks


def pymatgen_peaks(cif_path, wavelength=1.5406, two_theta_range=(10.0, 120.0)):
    from pymatgen.core import Structure
    from pymatgen.analysis.diffraction.xrd import XRDCalculator

    st = Structure.from_file(cif_path)
    pat = XRDCalculator(wavelength=wavelength).get_pattern(st, scaled=True, two_theta_range=two_theta_range)
    return pat, st


def pymatgen_form_factor(symbol, s):
    """pymatgen's own f(s) = Z - 41.78214 s^2 Σ a_i exp(-b_i s^2) (diagnostic use only)."""
    import json, os
    import pymatgen.analysis.diffraction as d
    from pymatgen.core import Element

    params = json.load(open(os.path.join(os.path.dirname(d.__file__), "atomic_scattering_params.json")))
    c = np.array(params[symbol])
    s2 = np.asarray(s, dtype=float) ** 2
    return Element(symbol).Z - 41.78214 * s2 * np.sum(c[:, 0] * np.exp(-np.multiply.outer(s2, c[:, 1])), axis=-1)


def compare(crystal, cif_path, wavelength=1.5406, two_theta_range=(10.0, 120.0), form_factor=None):
    kw = {} if form_factor is None else {"form_factor": form_factor}
    own = powder_peaks(crystal, wavelength, two_theta_range, **kw)
    pat, _ = pymatgen_peaks(cif_path, wavelength, two_theta_range)
    rows = []
    for p in own:
        k = int(np.argmin(np.abs(np.array(pat.x) - p.two_theta)))
        rows.append({
            "hkl": "+".join("".join(str(x) for x in f) for f in p.families),
            "d_own": p.d, "tt_own": p.two_theta, "tt_pmg": float(pat.x[k]),
            "dtt": p.two_theta - float(pat.x[k]),
            "I_own": p.intensity, "I_pmg": float(pat.y[k]), "dI": p.intensity - float(pat.y[k]),
            "mult": p.multiplicity,
        })
    return rows, own, pat
