"""最小 CIF 解析器 (minimal CIF reader) —— 只支持本项目需要的字段.

读取: 晶胞参数 _cell_length_*/_cell_angle_*, 对称操作
(_space_group_symop_operation_xyz 或 _symmetry_equiv_pos_as_xyz),
原子位置 _atom_site_label/_atom_site_fract_* (可选 _atom_site_type_symbol,
_atom_site_occupancy)。对不对称单元原子施加全部对称操作, 取模 1 并去重,
得到晶胞内全部原子 (full unit cell)。不依赖 pymatgen/ASE。
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np


@dataclass
class Crystal:
    name: str
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float
    species: list[str] = field(default_factory=list)      # one per atom in full cell
    frac_coords: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    occupancies: np.ndarray = field(default_factory=lambda: np.zeros(0))
    spacegroup: str = ""

    @property
    def cell(self):
        return (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)

    @property
    def n_atoms(self):
        return len(self.species)


def _num(s: str) -> float:
    """Parse CIF number, dropping standard uncertainty in parentheses: 5.43(1) -> 5.43."""
    return float(re.sub(r"\(.*?\)", "", s.strip().strip("'\"")))


def _parse_symop(op: str):
    """'1/4-x,3/4+z,y' -> (R 3x3 int, t 3 float) such that x' = R x + t."""
    op = op.strip().strip("'\"").replace(" ", "").lower()
    parts = op.split(",")
    if len(parts) != 3:
        raise ValueError(f"bad symop {op!r}")
    R = np.zeros((3, 3))
    t = np.zeros(3)
    for i, expr in enumerate(parts):
        for sign, body in re.findall(r"([+-]?)([^+-]+)", expr):
            s = -1.0 if sign == "-" else 1.0
            if body in "xyz":
                R[i, "xyz".index(body)] += s
            else:
                t[i] += s * float(Fraction(body))
    return R, t


def read_cif(path: str, name: str | None = None, dedupe_tol: float = 1e-4) -> Crystal:
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    scalars: dict[str, str] = {}
    loops: list[tuple[list[str], list[list[str]]]] = []
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln.startswith(";"):  # multi-line text field: skip
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                i += 1
            i += 1
            continue
        if ln.lower() == "loop_":
            i += 1
            tags = []
            while i < len(lines) and lines[i].strip().startswith("_"):
                tags.append(lines[i].strip().split()[0].lower())
                i += 1
            rows = []
            while i < len(lines):
                s = lines[i].strip()
                if not s or s.startswith("_") or s.lower() == "loop_" or s.startswith("data_") or s.startswith("#"):
                    break
                # symop lines may contain no spaces; atom lines are whitespace-separated
                rows.append(re.findall(r"'[^']*'|\"[^\"]*\"|\S+", s))
                i += 1
            loops.append((tags, rows))
            continue
        if ln.startswith("_"):
            toks = ln.split(None, 1)
            if len(toks) == 2:
                scalars[toks[0].lower()] = toks[1]
        i += 1

    a, b, c = (_num(scalars[f"_cell_length_{x}"]) for x in "abc")
    al, be, ga = (_num(scalars[f"_cell_angle_{x}"]) for x in ("alpha", "beta", "gamma"))
    sg = scalars.get("_symmetry_space_group_name_h-m", scalars.get("_space_group_name_h-m_alt", "")).strip("'\"")

    ops = None
    atoms = None
    for tags, rows in loops:
        for key in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"):
            if key in tags:
                j = tags.index(key)
                # a symop row tokenised by whitespace: rejoin everything after leading id columns
                ops = [_parse_symop("".join(r[j:]) if len(tags) == 1 else r[j]) for r in rows]
        if "_atom_site_fract_x" in tags:
            atoms = (tags, rows)
    if ops is None:
        ops = [(np.eye(3), np.zeros(3))]
    if atoms is None:
        raise ValueError("no _atom_site loop")

    tags, rows = atoms
    ix, iy, iz = (tags.index(f"_atom_site_fract_{q}") for q in "xyz")
    itype = tags.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in tags else tags.index("_atom_site_label")
    iocc = tags.index("_atom_site_occupancy") if "_atom_site_occupancy" in tags else None

    species, coords, occs = [], [], []
    for r in rows:
        el = re.match(r"[A-Z][a-z]?", r[itype].strip("'\"")).group(0)
        x = np.array([_num(r[ix]), _num(r[iy]), _num(r[iz])])
        occ = _num(r[iocc]) if iocc is not None else 1.0
        for R, t in ops:
            y = np.mod(R @ x + t, 1.0)
            y[np.isclose(y, 1.0, atol=dedupe_tol)] = 0.0
            dup = False
            for s2, c2 in zip(species, coords):
                if s2 == el:
                    dd = y - c2
                    dd -= np.round(dd)
                    if np.all(np.abs(dd) < dedupe_tol):
                        dup = True
                        break
            if not dup:
                species.append(el)
                coords.append(y)
                occs.append(occ)
    return Crystal(name or path, a, b, c, al, be, ga, species, np.array(coords), np.array(occs), sg)
