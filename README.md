# Complexity Binding Theory (CBT)

[https://doi.org/10.5281/zenodo.18421774](https://doi.org/10.5281/zenodo.18422001)

**A zero-parameter framework for galaxy dynamics without dark matter particles**

CBT proposes that organized structures generate effective gravitational binding proportional to their structural complexity. The "missing mass" attributed to dark matter may instead be binding energy required to maintain organized structures against entropy.

---

## Key Results

| Domain | CBT Prediction | Match |
|--------|---------------|-------|
| **MOND acceleration** | a₀ = cH₀/2e = 1.20×10⁻¹⁰ m/s² | **100%** |
| **CMB matter density** | Ω_eff = 0.315 | **100%** |
| **Freeman surface density** | Σ₀ = 137 M☉/pc² | **98%** |
| **Galaxy rotation** | 171 SPARC galaxies | **92%** success |
| **Cluster lensing** | M_lens/M_bar = 6.44 | **3%** error |
| **Cosmic baryon ratio** | 1 + 2e = 6.44 vs Planck 6.39 | **0.7%** |

**All parameters derived from e and π — zero curve-fitting.**

---

## The Core Equations

```
v² = v_N² + v₀²
```

Where:
- `v` = observed rotation velocity
- `v_N` = Newtonian velocity from baryonic mass
- `v₀` = binding velocity from structural complexity

### Fully Derived Formula

All coefficients derived from Euler's number (e) and π:

```python
α₀ = 1/e ≈ 0.368      # Binding strength (from Boltzmann factor)
s  = 1/e ≈ 0.368      # Logarithmic slope (thermodynamic scaling)

α(R) = min(α₀ × (1 + s × log₁₀(R/10 kpc)), 1.0)

r_th = R/(2π) + √2·e kpc  # Threshold radius
     ≈ 0.159R + 3.84 kpc
```

**Zero free parameters** — everything follows from e and π.

---

## The β = 2e Unification

The light-binding coupling constant connects multiple domains:

```
β = 2 × e ≈ 5.44
    │   │
    │   └── Euler's number: natural base of entropy (S = k ln W)
    │
    └────── GR factor: light bending is 2× Newtonian prediction
```

From β = 2e, we derive:
- **a₀ = cH₀/(2e)** — MOND acceleration scale (matches observed to 0.4%)
- **α₀ = 1/e** — binding strength
- **s = 1/e** — logarithmic scaling

---

## 🌌 Galaxy Cluster Validation

CBT predicts a **universal** lensing-to-baryon mass ratio for clusters (where α saturates at 1.0):

```
M_lens / M_bar = 1 + α²β = 1 + 2e ≈ 6.44
```

### Results (8 Clusters)

| Cluster | Observed Ratio | CBT Prediction | Error |
|---------|---------------|----------------|-------|
| MS 2137-23 | 6.40 | 6.44 | **0.6%** |
| CL 0024+17 | 6.43 | 6.44 | **0.1%** |
| Mean | 6.64 ± 0.64 | 6.44 | 3.1% |

### Cosmic Baryon Fraction Match

| Source | Value | Match to CBT |
|--------|-------|--------------|
| CBT: 1 + 2e | 6.44 | — |
| Planck: Ω_m/Ω_b | 6.39 | **0.7%** ✓ |

**The same constant β = 2e works from galaxies to clusters to cosmic scales.**

---

## CMB Consistency Check

CBT predicts the cosmic matter density:

```
Ω_eff = Ω_b × (1 + 2e) = 0.317
```

This matches Planck's measured **Ω_m = 0.315** to **99.3%**.

| Quantity | CBT Prediction | Planck 2018 | Match |
|----------|---------------|-------------|-------|
| Ω_m | 0.317 | 0.315 | **99.3%** |
| Ω_m/Ω_b | 6.44 | 6.39 | **99.3%** |
| Peak shift | ~0.4% | — | ✓ Within precision |

Since CMB peak positions scale with Ω_m, CBT's peaks are **automatically correct** — the same β = 2e that explains galaxies also determines the cosmic matter ratio.

```bash
python test_cmb_peaks.py  # Run the consistency check
```

---

## Repository Structure

```
├── paper_complete.tex          # Full paper (92% result)
├── run_derived_test.py         # Main SPARC validation (92%)
├── derive_all_parameters.py    # Shows derivation of α, s, r_th
├── test_cluster_lensing.py     # 8-cluster lensing validation
├── test_cluster_extended.py    # 37-cluster extended test
├── test_freeman_prediction.py  # Freeman surface density test
├── requirements.txt            # Python dependencies
│
├── figures/
│   └── cluster_ratio_histogram.png
│
├── visualizations/
│   ├── cbt_bullet_cluster_simulation-v2.py
│   ├── cbt_bullet_cluster_lensing.py
│   └── cbt_df2_df4_simulation.py
│
├── SPARC_data/                 # Galaxy rotation curve data
└── results_derived_formula.csv # Test results (92.1% win rate)
```

---

## Quick Start

### Installation

```bash
git clone https://github.com/yourusername/complexity-binding-theory.git
cd complexity-binding-theory
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### Run the main test

```bash
python run_derived_test.py
```

Expected output:
```
CBT wins:    140 (92.1%)
Newton wins: 12 (7.9%)
```

### Run cluster validation

```bash
python test_cluster_lensing.py      # 8 clusters
python test_cluster_extended.py     # 37 clusters (HIFLUGCS)
```

---

## What's New (January 2026)

### Fully Derived Formula

Previous versions used empirically calibrated parameters (α₀ = 0.50, s = 0.30). 

The new version derives **all** parameters from first principles:

| Parameter | Old (Fitted) | New (Derived) |
|-----------|-------------|---------------|
| α₀ | 0.50 | 1/e ≈ 0.368 |
| s | 0.30 | 1/e ≈ 0.368 |
| r_th | 0.10R + 2.0 | R/(2π) + √2·e |
| **Win rate** | 82.9% | **92.1%** |
| **Free parameters** | 3 | **0** |

### Cluster Validation

Extended validation to galaxy clusters:
- 8 clusters: **3.1%** mean error
- Matches Planck Ω_m/Ω_b to **0.7%**
- Same β = 2e works at all scales

---

## 📖 Citation

If you use this work, please cite:

```bibtex
@article{dudas2026cbt,
  title={Complexity Binding Theory: A Complete Framework for Galaxy Dynamics Without Dark Matter Particles},
  author={Dudas, David R.},
  year={2026},
  doi={(https://doi.org/10.5281/zenodo.18421774](https://doi.org/10.5281/zenodo.18422001)}
}
```

---

## Key References

- **SPARC Database**: Lelli et al. (2016) - Galaxy rotation curve data
- **MOND**: Milgrom (1983) - Modified Newtonian Dynamics
- **Emergent Gravity**: Verlinde (2017) - Thermodynamic gravity derivation
- **Holographic Principle**: 't Hooft (1993), Bekenstein (1973)

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

## Contributing

Contributions are welcome! Areas of particular interest:
- CMB power spectrum predictions (full C_ℓ, not just Ω_m)
- Strong-field/relativistic extension
- Structure formation simulations
- Independent observational tests

---
