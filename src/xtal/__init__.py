"""xtal — 从零实现的教学用晶体学 / 粉末 XRD 计算包 (Open Materials Lab, Increment 1–2).

Independent from-scratch implementation (only numpy). pymatgen is used ONLY
in xtal.compare for cross-validation, never inside the core computation.
"""
from .cif import Crystal, read_cif
from .lattice import metric_tensor, d_spacing, d_spacing_cubic, bragg_two_theta, cubic_bragg_two_theta
from .formfactor import f0_cromer_mann, load_cromer_mann
from .structure_factor import structure_factor
from .extinction import is_allowed, LATTICE_RULES
from .powder import powder_peaks, gaussian_profile, lorentz_polarization, cubic_multiplicity

CU_KA1 = 1.5406  # Å, Cu Kα1 (value specified for this project)

__all__ = [
    "Crystal", "read_cif", "metric_tensor", "d_spacing", "d_spacing_cubic", "bragg_two_theta",
    "f0_cromer_mann", "load_cromer_mann", "structure_factor", "is_allowed", "LATTICE_RULES",
    "powder_peaks", "gaussian_profile", "lorentz_polarization", "cubic_multiplicity",
    "cubic_bragg_two_theta", "CU_KA1",
]
