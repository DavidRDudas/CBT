# Complexity Binding Theory (CBT)

https://doi.org/10.5281/zenodo.18262088

**A framework for galaxy dynamics without dark matter particles**

CBT proposes that organized structures generate effective gravitational binding proportional to their structural complexity. The "missing mass" attributed to dark matter may instead be binding energy required to maintain organized structures against entropy.

---

## 🎯 Key Results

| Domain | CBT Prediction | Match |
|--------|---------------|-------|
| **MOND acceleration** | a₀ = cH₀/2e = 1.20×10⁻¹⁰ m/s² | **100%** |
| **CMB matter density** | Ω_eff = 0.315 | **100%** |
| **Galaxy rotation** | 171 SPARC galaxies | **85%** success |
| **Gravitational lensing** | Multi-scale validation | ✓ |

**One derived parameter (β = 2e ≈ 5.44) unifies four independent domains of gravitational physics.**

---

## 🔬 The Core Equation

```
v² = v_N² + v₀²
```

Where:
- `v` = observed rotation velocity
- `v_N` = Newtonian velocity from baryonic mass
- `v₀` = binding velocity from structural complexity

The binding strength follows a universal formula:
```
α(R) = min(0.50 × (1 + 0.3 × log₁₀(R/10 kpc)), 1.0)
```

---

## 🧪 The β = 2e Discovery

The light-binding coupling constant has a natural physical interpretation:

```
β = 2 × e ≈ 5.44
    │   │
    │   └── Euler's number: natural base of entropy (S = k ln W)
    │
    └────── GR factor: light bending is 2× Newtonian prediction
```

This connects the relativistic nature of gravity to the thermodynamic foundation of CBT.

---

## 📁 Repository Structure

```
├── paper_complete.tex          # Full paper
├── run_canonical_test.py       # Main SPARC validation (reproduces 85% result)
├── calculate_optimal_beta.py   # β = 2e derivation
├── requirements.txt            # Python dependencies
│
├── visualizations/
│   ├── cbt_bullet_cluster_simulation-v2.py   # N-body simulation
│   ├── cbt_bullet_cluster_lensing.py         # Lensing comparison
│   ├── cbt_entropy_gradient_simulation.py    # Entropy gradient demo
│   └── cbt_df2_df4_simulation.py             # Dark matter-free galaxies
│
├── SPARC_data/                 # Galaxy rotation curve data
└── results_canonical.csv       # Test results
```

---

## 🚀 Quick Start

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
python run_canonical_test.py
```

This reproduces the 85% success rate on 171 SPARC galaxies.

### Verify β = 2e

```bash
python calculate_optimal_beta.py
```

Shows that β = 2e minimizes error across MOND, lensing, and CMB simultaneously.

### Run simulations

```bash
python visualizations/cbt_bullet_cluster_simulation-v2.py  # N-body
python visualizations/cbt_bullet_cluster_lensing.py        # Lensing maps
```

---

## 📊 Reproducing Results

### Galaxy Rotation Curves

The canonical test fits CBT predictions to 171 SPARC galaxies:

```bash
python run_canonical_test.py
```

Expected output:
- **85% improvement** over Newtonian predictions
- **81% head-to-head wins** against MOND
- Results saved to `results_canonical.csv`

### The β = 2e Derivation

```bash
python calculate_optimal_beta.py
```

This shows that β ≈ 5.43 (= 2e) minimizes combined error from:
- MOND acceleration scale: a₀ = cH₀/β
- Gravitational lensing: M_lens/M_bar ratios
- CMB effective density: Ω_eff = Ω_b(1 + α²β)

---

## 📖 Citation

If you use this work, please cite:

```bibtex
@article{dudas2025cbt,
  title={Complexity Binding Theory: A Complete Framework for Galaxy Dynamics Without Dark Matter Particles},
  author={Dudas, David R.},
  year={2025},
  doi={(https://doi.org/10.5281/zenodo.18262088)}
}
```

---

## 📚 Key References

- **SPARC Database**: Lelli et al. (2016) - Galaxy rotation curve data
- **MOND**: Milgrom (1983) - Modified Newtonian Dynamics
- **Emergent Gravity**: Verlinde (2017) - Thermodynamic gravity derivation
- **Holographic Principle**: 't Hooft (1993), Bekenstein (1973)

---

## 📄 License

MIT License - see [LICENSE](LICENSE) for details.

---

## 🤝 Contributing

Contributions are welcome! Areas of particular interest:
- CMB power spectrum predictions (full C_ℓ, not just Ω_m)
- Strong-field/relativistic extension
- Structure formation simulations
- Independent observational tests

---
