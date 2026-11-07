"""晶格几何: 度规张量 (metric tensor), 晶面间距 d_hkl, Bragg 角.

一般情形 (any crystal system):
    G = [[a·a, a·b, a·c], ...]           (real-space metric tensor)
    1/d_hkl^2 = h^T G^{-1} h             (G^{-1} = reciprocal metric tensor)
立方 (cubic): d_hkl = a / sqrt(h^2 + k^2 + l^2)
Bragg 定律: λ = 2 d sinθ  ->  2θ = 2 arcsin(λ / 2d)
"""
from __future__ import annotations

import numpy as np


def metric_tensor(a, b, c, alpha, beta, gamma):
    al, be, ga = np.radians([alpha, beta, gamma])
    return np.array([
        [a * a, a * b * np.cos(ga), a * c * np.cos(be)],
        [a * b * np.cos(ga), b * b, b * c * np.cos(al)],
        [a * c * np.cos(be), b * c * np.cos(al), c * c],
    ])


def d_spacing(hkl, cell):
    """General d_hkl (Å) from the reciprocal metric tensor. cell=(a,b,c,α,β,γ)."""
    G_star = np.linalg.inv(metric_tensor(*cell))
    h = np.asarray(hkl, dtype=float)
    q2 = np.einsum("...i,ij,...j->...", h, G_star, h)
    return 1.0 / np.sqrt(q2)


def d_spacing_cubic(hkl, a):
    h = np.asarray(hkl, dtype=float)
    return a / np.sqrt(np.sum(h * h, axis=-1))


def bragg_two_theta(d, wavelength=1.5406):
    """2θ in degrees; NaN if λ/2d > 1 (reflection not accessible)."""
    s = wavelength / (2.0 * np.asarray(d, dtype=float))
    with np.errstate(invalid="ignore"):
        return np.degrees(2.0 * np.arcsin(np.where(s <= 1.0, s, np.nan)))


def cubic_bragg_two_theta(hkl, a, wavelength=1.5406):
    """Cubic Bragg 2θ (degrees) from lattice constant a (Å), Miller indices, and λ (Å).

    Convenience wrapper: d = a / sqrt(h²+k²+l²), then 2θ = 2 arcsin(λ / 2d).
    Returns NaN if the reflection is inaccessible (λ/2d > 1).
    """
    return bragg_two_theta(d_spacing_cubic(hkl, a), wavelength=wavelength)
