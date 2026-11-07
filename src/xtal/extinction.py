"""系统消光规则 (systematic-extinction / reflection conditions), 解析形式.

结构因子推导 (单原子 f, 晶胞内原子位置):
  BCC (Im-3m; 0,0,0 + ½,½,½):
      F = f[1 + (-1)^(h+k+l)]   -> 允许 h+k+l 偶数, F=2f
  FCC (Fm-3m; 0,0,0 + ½,½,0 + ½,0,½ + 0,½,½):
      F = f[1 + (-1)^(h+k) + (-1)^(h+l) + (-1)^(k+l)]
        -> h,k,l 全奇或全偶: F=4f; 混合奇偶: F=0
  Diamond (Fd-3m; FCC 格子 + 基元 (0,0,0),(¼,¼,¼)):
      F = F_FCC · [1 + i^(h+k+l)]
        -> 全奇: |F| = 4√2 f ; 全偶且 h+k+l=4n: |F| = 8f ;
           全偶且 h+k+l=4n+2: F = 0 (如 200, 222) ; 混合: F = 0
  注: 在本模型 (球形原子, 无成键电荷密度) 中 Si(222) 严格为零;
      真实 Si 晶体中可观察到很弱的 "forbidden" 222 反射 (源于成键电荷分布的非球形及
      热振动非谐性), 这是简单模型的已知局限。
"""
from __future__ import annotations

import numpy as np


def _parity(hkl):
    h = np.asarray(hkl, dtype=int)
    return h % 2


def allowed_bcc(hkl):
    return int(np.sum(hkl)) % 2 == 0


def allowed_fcc(hkl):
    p = _parity(hkl)
    return bool(np.all(p == p[0]))


def allowed_diamond(hkl):
    if not allowed_fcc(hkl):
        return False
    p = _parity(hkl)
    if p[0] == 1:  # all odd
        return True
    return int(np.sum(hkl)) % 4 == 0


def allowed_sc(hkl):
    return True


LATTICE_RULES = {"sc": allowed_sc, "bcc": allowed_bcc, "fcc": allowed_fcc, "diamond": allowed_diamond}


def is_allowed(hkl, kind: str) -> bool:
    return LATTICE_RULES[kind](tuple(int(x) for x in hkl))


def analytic_F_over_f(hkl, kind: str) -> complex:
    """|F|/f predicted analytically (geometric part only)."""
    h, k, l = (int(x) for x in hkl)
    if kind == "bcc":
        return 1 + (-1) ** (h + k + l)
    fcc = 1 + (-1) ** (h + k) + (-1) ** (h + l) + (-1) ** (k + l)
    if kind == "fcc":
        return fcc
    if kind == "diamond":
        return fcc * (1 + 1j ** (h + k + l))
    return 1
