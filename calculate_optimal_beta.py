#!/usr/bin/env python3
"""
Calculate the optimal value of β from multiple independent constraints.

Constraints:
1. MOND acceleration: a₀ = cH₀/β should match ~1.2×10⁻¹⁰ m/s²
2. Lensing: M_lens/M_bar ratios for galaxies and clusters
3. CMB: Ω_eff = Ω_b(1 + α²β) should match Planck's Ω_m = 0.315
"""

import numpy as np

# Physical constants
c = 2.998e8  # m/s
H0_planck = 67.4  # km/s/Mpc (Planck 2018)
H0_si = H0_planck * 1000 / (3.086e22)  # Convert to s⁻¹
print(f"H₀ = {H0_si:.3e} s⁻¹")
print(f"cH₀ = {c * H0_si:.3e} m/s²")

# =============================================================================
# CONSTRAINT 1: MOND ACCELERATION SCALE
# =============================================================================
a0_MOND = 1.2e-10  # m/s² (empirical, from galaxy fitting)
cH0 = c * H0_si

print("\n" + "="*60)
print("CONSTRAINT 1: MOND ACCELERATION SCALE")
print("="*60)
print(f"Target: a₀_MOND = {a0_MOND:.2e} m/s²")

for beta in [5, 6, 7, 8, 9]:
    a0_pred = cH0 / beta
    match = a0_pred / a0_MOND * 100
    print(f"β = {beta}: a₀ = cH₀/{beta} = {a0_pred:.2e} m/s² ({match:.1f}% of MOND)")

beta_optimal_mond = cH0 / a0_MOND
print(f"\nOptimal β for MOND: {beta_optimal_mond:.2f}")

# =============================================================================
# CONSTRAINT 2: GRAVITATIONAL LENSING
# =============================================================================
print("\n" + "="*60)
print("CONSTRAINT 2: GRAVITATIONAL LENSING")
print("="*60)

# Typical observations (from various surveys)
# Galaxies: dynamical mass / baryonic mass ≈ 2-5
# Clusters: lensing mass / baryonic mass ≈ 5-8

# For a Milky Way type galaxy: α ≈ 0.53, observed M_dyn/M_bar ≈ 3
alpha_galaxy = 0.53
observed_ratio_galaxy = 3.0

# For a cluster: α = 1.0 (saturated), observed M_lens/M_bar ≈ 6
alpha_cluster = 1.0
observed_ratio_cluster = 6.0

print(f"\nGalaxy (α={alpha_galaxy}, observed M/M_bar={observed_ratio_galaxy}):")
for beta in [5, 6, 7, 8, 9]:
    pred = 1 + alpha_galaxy**2 * beta
    error = abs(pred - observed_ratio_galaxy) / observed_ratio_galaxy * 100
    print(f"  β = {beta}: M_lens/M_bar = {pred:.2f} (error: {error:.1f}%)")

print(f"\nCluster (α={alpha_cluster}, observed M/M_bar={observed_ratio_cluster}):")
for beta in [5, 6, 7, 8, 9]:
    pred = 1 + alpha_cluster**2 * beta
    error = abs(pred - observed_ratio_cluster) / observed_ratio_cluster * 100
    print(f"  β = {beta}: M_lens/M_bar = {pred:.2f} (error: {error:.1f}%)")

# Optimal β for lensing (average of both)
beta_optimal_galaxy = (observed_ratio_galaxy - 1) / alpha_galaxy**2
beta_optimal_cluster = (observed_ratio_cluster - 1) / alpha_cluster**2
print(f"\nOptimal β for galaxy lensing: {beta_optimal_galaxy:.2f}")
print(f"Optimal β for cluster lensing: {beta_optimal_cluster:.2f}")

# =============================================================================
# CONSTRAINT 3: CMB
# =============================================================================
print("\n" + "="*60)
print("CONSTRAINT 3: CMB EFFECTIVE MATTER DENSITY")
print("="*60)

Omega_b = 0.049  # Planck baryon density
Omega_m_planck = 0.315  # Planck total matter density
alpha_CMB = 1.0  # Saturated at sound horizon scale

print(f"Target: Ω_m = {Omega_m_planck}")
print(f"Formula: Ω_eff = Ω_b × (1 + α²β) = {Omega_b} × (1 + β)")

for beta in [5, 6, 7, 8, 9]:
    Omega_eff = Omega_b * (1 + alpha_CMB**2 * beta)
    match = Omega_eff / Omega_m_planck * 100
    print(f"  β = {beta}: Ω_eff = {Omega_eff:.3f} ({match:.1f}% of Planck)")

beta_optimal_CMB = (Omega_m_planck / Omega_b - 1) / alpha_CMB**2
print(f"\nOptimal β for CMB: {beta_optimal_CMB:.2f}")

# =============================================================================
# COMBINED OPTIMAL β
# =============================================================================
print("\n" + "="*60)
print("COMBINED ANALYSIS")
print("="*60)

# Calculate chi-squared-like metric for each β
print("\nTotal error for each β (sum of squared fractional errors):")
for beta in [4, 5, 5.5, 6, 6.5, 7, 8, 9]:
    # MOND error
    a0_pred = cH0 / beta
    err_mond = ((a0_pred - a0_MOND) / a0_MOND)**2
    
    # Lensing error (cluster, more reliable)
    lensing_pred = 1 + alpha_cluster**2 * beta
    err_lensing = ((lensing_pred - observed_ratio_cluster) / observed_ratio_cluster)**2
    
    # CMB error
    Omega_eff = Omega_b * (1 + alpha_CMB**2 * beta)
    err_CMB = ((Omega_eff - Omega_m_planck) / Omega_m_planck)**2
    
    total_err = err_mond + err_lensing + err_CMB
    print(f"  β = {beta:4.1f}: MOND={err_mond:.3f}, Lens={err_lensing:.3f}, CMB={err_CMB:.3f} → Total={total_err:.3f}")

# Optimize
from scipy.optimize import minimize_scalar

def total_error(beta):
    a0_pred = cH0 / beta
    err_mond = ((a0_pred - a0_MOND) / a0_MOND)**2
    lensing_pred = 1 + alpha_cluster**2 * beta
    err_lensing = ((lensing_pred - observed_ratio_cluster) / observed_ratio_cluster)**2
    Omega_eff = Omega_b * (1 + alpha_CMB**2 * beta)
    err_CMB = ((Omega_eff - Omega_m_planck) / Omega_m_planck)**2
    return err_mond + err_lensing + err_CMB

result = minimize_scalar(total_error, bounds=(3, 12), method='bounded')
print(f"\n*** OPTIMAL β (minimizing total error): {result.x:.2f} ***")
print(f"    (Integer approximation: {round(result.x)})")
