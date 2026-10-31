"""原子散射因子 f0 (atomic form factor), Cromer–Mann 9 参数形式:

    f0(s) = Σ_{i=1..4} a_i exp(-b_i s^2) + c ,   s = sinθ/λ  (Å^-1)

系数来源 (coefficient source):
  International Tables for Crystallography / for X-ray Crystallography, Vol. C
  (ed. A. J. C. Wilson, Kluwer, 1995), pp. 500–502 —— 通过 ESRF DABAX 库文件
  f0_InterTables.dat 获取 (https://ftp.esrf.fr/pub/scisoft/xop2.3/DabaxFiles/f0_InterTables.dat),
  副本保存在 data/dabax_f0_InterTables.dat。参数化有效范围 0 ≤ s < 2.0 Å^-1。
不包含反常色散 f', f'' (anomalous dispersion), 也不包含 Debye–Waller 因子。
"""
from __future__ import annotations

import os
import re
from functools import lru_cache

import numpy as np

_DEFAULT = os.path.join(os.path.dirname(__file__), "..", "..", "data", "dabax_f0_InterTables.dat")


@lru_cache(maxsize=None)
def load_cromer_mann(path: str = _DEFAULT) -> dict[str, tuple[np.ndarray, np.ndarray, float]]:
    """Return {symbol: (a[4], b[4], c)} for neutral atoms in a DABAX f0 file."""
    out = {}
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    for i, ln in enumerate(lines):
        m = re.match(r"#S\s+\d+\s+(\S+)\s*$", ln)
        if not m:
            continue
        sym = m.group(1)
        if not re.fullmatch(r"[A-Z][a-z]?", sym):
            continue  # skip ions like Fe+2, Si4+, 'Si.'
        j = i + 1
        while lines[j].startswith("#"):
            j += 1
        v = [float(x) for x in lines[j].split()]
        a1, a2, a3, a4, c, b1, b2, b3, b4 = v  # DABAX column order
        out[sym] = (np.array([a1, a2, a3, a4]), np.array([b1, b2, b3, b4]), c)
    return out


def f0_cromer_mann(symbol: str, s, path: str = _DEFAULT):
    a, b, c = load_cromer_mann(path)[symbol]
    s2 = np.asarray(s, dtype=float) ** 2
    return np.sum(a * np.exp(-np.multiply.outer(s2, b)), axis=-1) + c
