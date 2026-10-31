"""结构因子 (structure factor):

    F(hkl) = Σ_j occ_j · f_j(s) · exp[2πi (h x_j + k y_j + l z_j)]

Debye–Waller (温度) 因子默认关闭 (B = 0); 可选各向同性 B (Å^2):
    f_j -> f_j · exp(-B_j s^2)
"""
from __future__ import annotations

import numpy as np

from .formfactor import f0_cromer_mann


def structure_factor(crystal, hkl, s, form_factor=f0_cromer_mann, B: dict | None = None):
    """hkl: (N,3) ints; s: (N,) sinθ/λ. Returns complex F (N,)."""
    hkl = np.atleast_2d(np.asarray(hkl, dtype=float))
    s = np.atleast_1d(np.asarray(s, dtype=float))
    F = np.zeros(len(hkl), dtype=complex)
    for el in sorted(set(crystal.species)):
        mask = np.array([sp == el for sp in crystal.species])
        xyz = crystal.frac_coords[mask]
        occ = crystal.occupancies[mask]
        f = form_factor(el, s)
        if B and el in B:
            f = f * np.exp(-B[el] * s ** 2)
        phase = np.exp(2j * np.pi * hkl @ xyz.T)  # (N, n_el)
        F += f * (phase @ occ)
    return F
