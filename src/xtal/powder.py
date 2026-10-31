"""粉末衍射花样 (powder diffraction pattern), 教学模拟.

对每个 (hkl), 在极限球 |g| ≤ 2/λ 内枚举所有倒易格点, 按 d 值 (即 2θ) 合并为一个峰:
    I_peak ∝ Σ_{hkl ∈ peak} |F(hkl)|^2 · LP(θ)
求和遍历全部等价 hkl, 等价于 "多重性 p × |F|^2" (multiplicity)。
Lorentz–polarization 因子 (无单色器, 非偏振入射):
    LP(θ) = (1 + cos^2 2θ) / (sin^2 θ · cos θ)
(省略常数因子 1/2 或 1/8, 因为只比较相对强度。)
未包含: Debye–Waller 因子 (默认 B=0), 吸收, 择优取向, 仪器函数, Kα2, 反常色散.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

from .formfactor import f0_cromer_mann
from .lattice import bragg_two_theta, metric_tensor
from .structure_factor import structure_factor


def lorentz_polarization(two_theta_deg):
    th = np.radians(np.asarray(two_theta_deg) / 2.0)
    return (1 + np.cos(2 * th) ** 2) / (np.sin(th) ** 2 * np.cos(th))


def cubic_multiplicity(hkl) -> int:
    """Number of distinct (h,k,l) related by the m-3m point group (cubic Laue class)."""
    h = tuple(abs(int(x)) for x in hkl)
    perms = set(itertools.permutations(h))
    out = set()
    for p in perms:
        for signs in itertools.product((1, -1), repeat=3):
            out.add(tuple(s * x for s, x in zip(signs, p)))
    return len(out)


@dataclass
class Peak:
    two_theta: float
    d: float
    intensity: float          # relative, max = 100
    raw_intensity: float      # Σ|F|^2·LP (unnormalised)
    families: list            # canonical hkl families merged into this peak (|h|>=|k|>=|l|)
    multiplicity: int         # total number of contributing hkl with non-zero F
    F_abs: float              # |F| of one representative reflection


def powder_peaks(crystal, wavelength=1.5406, two_theta_range=(10.0, 120.0), form_factor=f0_cromer_mann,
                 B=None, merge_tol=1e-5, min_rel_intensity=1e-3, zero_tol=1e-8):
    G_star = np.linalg.inv(metric_tensor(*crystal.cell))
    q_max = 2.0 / wavelength  # limiting sphere radius (1/d)
    # bound on |h| from reciprocal lattice: |h_i| <= q_max * |a_i| (a_i real lattice vector length)
    lens = np.array(crystal.cell[:3])
    hmax = np.ceil(q_max * lens).astype(int)
    rng = [np.arange(-m, m + 1) for m in hmax]
    H = np.array(np.meshgrid(*rng, indexing="ij")).reshape(3, -1).T
    H = H[np.any(H != 0, axis=1)]
    q2 = np.einsum("ni,ij,nj->n", H, G_star, H)
    d = 1 / np.sqrt(q2)
    tt = bragg_two_theta(d, wavelength)
    keep = np.isfinite(tt) & (tt >= two_theta_range[0]) & (tt <= two_theta_range[1])
    H, d, tt = H[keep], d[keep], tt[keep]
    s = 1.0 / (2.0 * d)
    F = structure_factor(crystal, H, s, form_factor=form_factor, B=B)
    F2 = np.abs(F) ** 2
    lp = lorentz_polarization(tt)

    order = np.argsort(tt)
    H, d, tt, F, F2, lp = H[order], d[order], tt[order], F[order], F2[order], lp[order]
    peaks = []
    i = 0
    scale = max(F2.max(), 1.0)
    while i < len(tt):
        j = i
        while j + 1 < len(tt) and abs(tt[j + 1] - tt[i]) < merge_tol:
            j += 1
        sel = slice(i, j + 1)
        nz = F2[sel] > zero_tol * scale
        if nz.any():
            fams = sorted({tuple(sorted(np.abs(h), reverse=True)) for h in H[sel][nz]}, reverse=True)
            raw = float(np.sum(F2[sel][nz] * lp[sel][nz]))
            peaks.append(Peak(float(tt[i]), float(d[i]), 0.0, raw, [tuple(int(x) for x in f) for f in fams],
                              int(nz.sum()), float(np.abs(F[sel][nz][0]))))
        i = j + 1
    Imax = max(p.raw_intensity for p in peaks)
    for p in peaks:
        p.intensity = 100.0 * p.raw_intensity / Imax
    return [p for p in peaks if p.intensity >= min_rel_intensity]


def gaussian_profile(peaks, two_theta_grid, fwhm=0.15):
    """Sum of Gaussians (fixed FWHM in degrees 2θ) — purely illustrative broadening."""
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    x = np.asarray(two_theta_grid)[:, None]
    c = np.array([p.two_theta for p in peaks])[None, :]
    h = np.array([p.intensity for p in peaks])[None, :]
    y = np.sum(h * np.exp(-0.5 * ((x - c) / sigma) ** 2), axis=1)
    return y
