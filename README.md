# Open Materials Lab — Increment 1–2：粉末 XRD 模拟 + NIST SRM 640 基准

> **教学模拟 (teaching simulation)，不是实验测量 (not experimental measurement)。**
> 本仓库是对教科书级晶体学计算的**可复现实现与交叉验证**，不包含任何原创性科研结论。

## 目的 (Purpose)

为本科生材料计算训练搭建第一个最小可运行、可验证的案例：
从公开数据库 **COD (Crystallography Open Database)** 下载 Si、Al、α-Fe 的实测晶体结构 (CIF)，
用**从零编写的 Python 代码** (`src/xtal/`, 仅依赖 numpy) 计算：

- 晶面间距 *d*-spacing（立方公式 + 通用度规张量 metric tensor）
- Bragg 衍射角 2θ（Cu Kα1, λ = 1.5406 Å）
- 原子散射因子 atomic form factor f₀(s)（Cromer–Mann 参数，International Tables Vol. C）
- 结构因子 structure factor F(hkl) 与 FCC / BCC / 金刚石结构的系统消光 (systematic extinctions)
- 多重性 multiplicity、Lorentz–polarization 因子、粉末花样 (stick + 可选 Gaussian 展宽)

并用成熟工具 **pymatgen `XRDCalculator`** 交叉验证峰位和相对强度，用 **ASE** 检查 CIF 读取。

## 目录结构

```
open-materials-lab/
├── README.md
├── requirements.txt          # 直接依赖 (精确版本)
├── requirements-lock.txt     # pip freeze 全量锁定
├── pyproject.toml            # pytest 配置
├── data/
│   ├── SOURCES.md            # 数据来源、COD ID、URL、下载日期、SHA-256
│   ├── Si_COD_9011998.cif
│   ├── Al_COD_9012002.cif
│   ├── Fe_COD_9008536.cif
│   ├── dabax_f0_InterTables.dat   # Cromer–Mann 系数 (Int. Tables Vol. C, via ESRF DABAX)
│   └── nist/                 # NIST SRM 640e/640f 证书 PDF + JSON/CSV 转录 (Increment 2)
├── src/xtal/
│   ├── cif.py                # 最小 CIF 解析 + 对称操作展开
│   ├── lattice.py            # 度规张量, d_hkl, Bragg 角, cubic_bragg_two_theta
│   ├── formfactor.py         # f0(s) Cromer–Mann
│   ├── structure_factor.py   # F(hkl)
│   ├── extinction.py         # 解析消光规则与推导
│   ├── powder.py             # 多重性, LP 因子, 粉末花样, Gaussian 展宽
│   └── compare.py            # 与 pymatgen 对比 (唯一导入 pymatgen 的模块)
├── tests/
│   ├── test_xtal.py          # Increment 1 pytest
│   └── test_nist_srm640.py   # Increment 2 NIST SRM 640 基准测试
├── scripts/
│   ├── run_validation.py     # 一键生成验证表与图
│   ├── run_nist_srm640_benchmark.py  # NIST SRM 640 晶格/峰位基准 (Increment 2)
│   └── build_notebook.py     # 生成 notebook
├── docs/PROGRESS.md          # 增量进度与下一步
├── notebooks/01_xrd_si_al_fe.ipynb (+ .html)   # 已用 nbconvert 完整执行
├── figures/{Si,Al,Fe}_xrd_own_vs_pymatgen.png
├── figures/Si_nist_srm640e_peak_positions.png
└── results/validation_*.csv, validation_summary.md, nist_srm640_*.csv/.md
```

## 如何复现 (How to reproduce)

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m pytest -v                         # 运行全部测试
python scripts/run_validation.py            # 生成 results/ 与 figures/ (Increment 1)
python scripts/run_nist_srm640_benchmark.py # NIST SRM 640 基准 (Increment 2)
cd notebooks && jupyter nbconvert --to notebook --execute --inplace 01_xrd_si_al_fe.ipynb
```

开发环境：Python 3.13.5, Linux x86_64。

## 物理模型 (Model)

| 量 | 公式 |
|---|---|
| d-spacing | 1/d² = hᵀ G⁻¹ h；立方 d = a/√(h²+k²+l²) |
| Bragg | λ = 2d sinθ |
| form factor | f₀(s) = Σ aᵢ exp(−bᵢ s²) + c, s = sinθ/λ |
| structure factor | F = Σⱼ fⱼ exp[2πi(hxⱼ+kyⱼ+lzⱼ)] |
| powder intensity | I ∝ Σ_equiv \|F\|² · (1+cos²2θ)/(sin²θ cosθ) |

系统消光（单原子基元，解析推导见 `src/xtal/extinction.py`）：
- **BCC** (α-Fe)：h+k+l 偶数才出现 → 无 (100)、(111)
- **FCC** (Al)：h,k,l 全奇或全偶 → 无 (100)、(110)
- **Diamond** (Si)：全奇，或全偶且 h+k+l = 4n → 无 (200)、(222)

## 结果摘要 (Results)

2θ 范围 10–120°，Cu Kα1。完整表格见 `results/validation_summary.md`。

| 材料 | 峰数 own / pymatgen | max \|Δ2θ\| | max \|ΔI\| (0–100 标度) | 换用相同 f(s) 后 max \|ΔI\| |
|---|---|---|---|---|
| Si | 9 / 9 | 2.8×10⁻¹⁴ ° | 0.45 (220) | 2.8×10⁻¹⁴ |
| Al | 8 / 8 | 5.7×10⁻¹⁴ ° | 0.78 (311) | 2.1×10⁻¹⁴ |
| α-Fe | 5 / 5 | 4.3×10⁻¹⁴ ° | 0.46 (211) | 1.8×10⁻¹⁴ |

前几个峰：

| 材料 | hkl | 2θ own (°) | 2θ pmg (°) | I own | I pmg |
|---|---|---|---|---|---|
| Si | 111 | 28.4424 | 28.4424 | 100.00 | 100.00 |
| Si | 220 | 47.3028 | 47.3028 | 67.11 | 66.66 |
| Si | 311 | 56.1225 | 56.1225 | 39.94 | 39.59 |
| Si | 400 | 69.1301 | 69.1301 | 10.70 | 10.71 |
| Al | 111 | 38.4724 | 38.4724 | 100.00 | 100.00 |
| Al | 200 | 44.7209 | 44.7209 | 48.37 | 48.14 |
| Al | 220 | 65.0968 | 65.0968 | 30.70 | 31.05 |
| Al | 311 | 78.2294 | 78.2294 | 35.54 | 36.31 |
| Fe | 110 | 44.6716 | 44.6716 | 100.00 | 100.00 |
| Fe | 200 | 65.0204 | 65.0204 | 14.65 | 14.69 |
| Fe | 211 | 82.3314 | 82.3314 | 28.02 | 28.49 |
| Fe | 220 | 98.9405 | 98.9405 | 9.61 | 9.75 |

**解读**：
1. 峰位一致到浮点舍入：两者都由同一 CIF 晶格常数和 Bragg 定律决定，这验证了 d/2θ 和 hkl 枚举，但不是对"实验峰位"的验证。
2. 强度差异 < 1（满分 100），来源是**原子散射因子表不同**（本项目用 Int. Tables Vol. C 的 Cromer–Mann 参数；pymatgen 用自带的另一套参数与 Mott–Bethe 形式公式）。把 pymatgen 的 f(s) 代入本代码后差异降至 ~10⁻¹⁴，说明 F(hkl)、多重性、LP 因子实现一致。
3. 消光：对 0 ≤ h,k,l ≤ 5 的全部 hkl，数值 |F| 与解析公式一致，且非零 ⇔ 解析规则判定"允许"。


## Increment 2 — NIST SRM 640 Si 晶格 / 峰位基准 (lattice & peak-position benchmark)

> **教学模拟 + 文献/COD 与 NIST 认证标准对比。** 本项目**没有**测量 SRM 粉末；证书 Table 1 为 NIST *information values* 的转录。不声称创新性。

对照 NIST SRM **640e** (证书 2015-10-29) 与 **640f** (证书 2019-02-26) 的认证晶格常数，以及证书 Table 1（Cu Kα λ = 0.15405929 nm）峰位信息值。证书 PDF、SHA-256、机读 JSON/CSV 见 `data/nist/` 与 `data/SOURCES.md`。完整表格见 `results/nist_srm640_benchmark.md`。

| 量 | 数值 |
|---|---|
| COD Si a (9011998, 298.15 K) | 5.430941 Å |
| NIST 640e certified a (22.5 °C) | 5.431179 Å ± 0.00008 Å |
| NIST 640f certified a (22.5 °C) | 5.431144 Å ± 0.00008 Å |
| Δa (COD − 640e) | −0.000238 Å (**−43.82 ppm**) |
| Δa (COD − 640f) | −0.000203 Å (**−37.38 ppm**) |
| max \|Δ2θ\| 用 640e 认证 a + 证书 λ 重算 Table 1 | < 0.001° (~0.00048°) |
| max \|Δ2θ\| 用 640f 认证 a + 证书 λ 重算 Table 1 | < 0.001° (~0.00050°) |

解读：
1. 用 NIST 认证 a 与证书波长重算 Bragg 2θ，与证书 Table 1 在 0.001° 量级内一致 —— 验证立方 Bragg 实现与证书转录，**不是**对 SRM 粉末的实验测量。
2. COD a 与 NIST 认证 a 相差约 40 ppm 量级，且温度略有不同（298.15 K vs 22.5 °C）；未做热膨胀修正。
3. 另用 COD a + 本项目 λ=1.5406 Å 对比 640e Table 1，仅作峰位灵敏度的教学示意（混合了不同 a 与略不同的 λ）。

进度与下一步见 [`docs/PROGRESS.md`](docs/PROGRESS.md)。下一步拟议：**各向同性 Debye–Waller B 因子**（文献引用的 B 值）。

## 局限 (Limitations)

- **教学模拟**：无 Debye–Waller (B = 0)、无吸收、无择优取向、无 Kα2、无仪器函数、无反常色散 f′/f″；Gaussian 展宽 FWHM = 0.15° 只是示意。高角峰的相对强度因此系统性偏高。
- 球形原子模型中 Si(222) 严格为零；真实 Si 晶体中可观察到很弱的 (222) 反射（源于成键电荷分布与非谐热振动），本模型不能描述。
- 峰位与 pymatgen 一致只证明代码一致；与 NIST SRM 640 认证晶格常数 / Table 1 峰位的对比见 **Increment 2**（仍非本实验室对 SRM 粉末的测量）。
- 自写 CIF 解析器只覆盖本项目所需字段，不是通用 CIF 解析器。
- 只验证了立方晶系的数值结果；通用度规张量公式仅用一个六方算例做了单元测试。

## 来源 (Sources)

详见 [`data/SOURCES.md`](data/SOURCES.md)：
- COD 9011998 (Si; Hom, Kiszenick & Post 1975)、9012002 (Al; Straumanis 1949)、9008536 (α-Fe; Wyckoff 1963)，https://www.crystallography.net/cod/
- Cromer–Mann 系数：International Tables for X-ray Crystallography Vol. C (1995)，经 ESRF DABAX `f0_InterTables.dat`
- pymatgen：Ong et al., Comput. Mater. Sci. 68, 314 (2013)
- NIST SRM 640e / 640f 证书（下载 2026-11-07）：https://tsapps.nist.gov/srmext/certificates/archives/640e.pdf 、…/640f.pdf

## 许可与诚信声明

代码为教学复现；不伪造数据，不声称创新性。CIF 文件遵循 COD 的使用条款并需引用原始文献。
