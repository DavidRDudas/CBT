# Complexity Binding Theory

A phenomenological framework for galaxy dynamics without dark matter particles.

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run the canonical SPARC evaluation
python run_canonical_test.py

# Or with BTF-predicted Vmax (stricter test)
python run_canonical_test.py --use-btf-vmax
```

## Parameter Accounting

| Model | Free Parameters |
|-------|-----------------|
| Newton | 1 (M/L scale) |
| CBT | 1 (M/L scale) + universal α(R) formula |

The α(R) formula is **universal and fixed**:
```
α(R) = 0.50 × (1 + 0.3 × log₁₀(R/10 kpc))
```

There is **no per-galaxy tuning** of α.

## Core Equation

```
v² = v_N² + v₀²
```

Where:
- v = observed rotation velocity
- v_N = Newtonian prediction from baryonic mass
- v₀ = binding velocity from complexity

## Repository Structure

```
├── run_canonical_test.py    # Official evaluation script
├── SPARC_data/              # Galaxy rotation curve data
├── paper_complete.tex       # The academic paper
├── requirements.txt         # Python dependencies
└── results_canonical.csv    # Output from evaluation
```

## Key Results

- **81%** of SPARC galaxies fit better with CBT than Newton
- **81%** head-to-head wins against MOND
- **104%** match on gravitational lensing masses
- **5 of 6** unique predictions supported by literature

## Data

SPARC data from Lelli et al. (2016): http://astroweb.cwru.edu/SPARC/

## Citation

If you use this code, please cite the paper:
```
Dudas, D. R. (2026). Complexity Binding Theory: A Complete Framework 
for Galaxy Dynamics Without Dark Matter Particles. arXiv:TBD
```

## License

MIT License - see LICENSE file.
