# 数据来源 (Data sources)

下载日期 (download date): **2026-10-31** (Asia/Shanghai, UTC+8)，通过 `curl` 从 COD 公开接口获取，无需 API key。
COD 文件说明：数据由 American Mineralogist Crystal Structure Database (AMCSD) 提供给 COD；使用时需引用原始文献。

## 晶体结构 (CIF)

| 文件 | 材料 | COD ID | URL | 空间群 | a (Å) | 温度 | 原始文献 | SHA-256 |
|---|---|---|---|---|---|---|---|---|
| `Si_COD_9011998.cif` | Si, diamond cubic | 9011998 (AMCSD 0012843) | https://www.crystallography.net/cod/9011998.cif | F d -3 m (origin choice 1), No. 227 | 5.430941 | 298.15 K | Hom, T.; Kiszenick, W.; Post, B. *J. Appl. Cryst.* **8**, 457–458 (1975). doi:10.1107/S0021889875010965 | a3dfc7a2fcd58df8f6f68215b182c1920b0d2e9a9cc15f285b73333d051d5cb6 |
| `Al_COD_9012002.cif` | Al, FCC | 9012002 (AMCSD 0012872) | https://www.crystallography.net/cod/9012002.cif | F m -3 m, No. 225 | 4.04961 | 298.15 K | Straumanis, M. E. *J. Appl. Phys.* **20**, 726–734 (1949). doi:10.1063/1.1698520 | d75379b3181af558a4a473cbeab3a6065ac83add6933301483a3f23a0aee1879 |
| `Fe_COD_9008536.cif` | α-Fe, BCC | 9008536 (AMCSD 0011214) | https://www.crystallography.net/cod/9008536.cif | I m -3 m, No. 229 | 2.8665 | 298 K | Wyckoff, R. W. G. *Crystal Structures*, 2nd ed., Vol. 1, pp. 7–83, Interscience, New York (1963) | 58d627740a3a7eae35f96e0a0f32be51c428792fb3b3d85f2215a9c082ad6591 |

备注：
- Al 条目 (9012002) 的 COD 维护者注释说明：原文晶胞常数为 kX 单位，COD 已乘以 1.00202 换算为 Å。
- 选择理由：均为室温 (≈298 K)、常压、非理论计算的实验测定条目；COD 中同一元素还有高压/高温条目，本项目未采用。
- 这些 CIF 只提供晶胞和原子位置，不含原子位移参数 (ADP / B 因子)，因此本项目 Debye–Waller 因子关闭。

## 原子散射因子系数 (atomic form factor coefficients)

- 文件：`dabax_f0_InterTables.dat`，来自 ESRF DABAX 库：https://ftp.esrf.fr/pub/scisoft/xop2.3/DabaxFiles/f0_InterTables.dat （下载日期 2026-10-31）
- 文件内注明的原始来源：*International Tables for X-ray Crystallography*, Vol. C (1995), ed. A. J. C. Wilson, Kluwer, pp. 500–502。
- 形式：f0(s) = Σ_{i=1..4} a_i exp(−b_i s²) + c，s = sinθ/λ，有效范围 0 ≤ s < 2.0 Å⁻¹。
- 本项目用到的中性原子系数 (a1 a2 a3 a4 c b1 b2 b3 b4)：
  - Al: 6.420200 1.900200 1.593600 1.964600 1.115100 3.038700 0.7426000 31.54720 85.08860
  - Si: 6.291500 3.035300 1.989100 1.541000 1.140700 2.438600 32.33370 0.6785000 81.69370
  - Fe: 11.76950 7.357300 3.522200 2.304500 1.036900 4.761100 0.3072000 15.35350 76.88050

## 交叉验证工具

- pymatgen 2026.10.24 `pymatgen.analysis.diffraction.xrd.XRDCalculator`（Ong et al., *Comput. Mater. Sci.* 68, 314–319 (2013), doi:10.1016/j.commatsci.2012.10.028）。
  其原子散射因子使用包内 `atomic_scattering_params.json`，公式 f(s) = Z − 41.78214·s²·Σ a_i exp(−b_i s²)（源码 docstring 称形式依据 De Graef & McHenry, *Structure of Materials* 第 11–12 章；源码未注明系数表的具体出处）。
- ASE 3.29.0：仅用于独立检查 CIF 读取的原子数与晶胞参数。
