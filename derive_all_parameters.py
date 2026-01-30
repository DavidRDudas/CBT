#!/usr/bin/env python3
"""
DERIVING ALL CBT PARAMETERS FROM FIRST PRINCIPLES
===================================================

Goal: Derive α₀, s, and r_th from β = 2e and physical principles.

Currently empirical:
- α₀ = 0.50
- s = 0.30  
- r_th = 0.10R + 2.0 kpc

Target: Express ALL in terms of e, π, or other fundamental constants.

Author: David R. Dudas
Date: January 2026
"""

import numpy as np
from scipy.optimize import curve_fit
import os
import glob
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# FUNDAMENTAL CONSTANTS
# =============================================================================

e = np.e
pi = np.pi
BETA = 2 * e

# Physical constants
c = 2.998e5  # km/s
H0 = 70  # km/s/Mpc
a0_SI = c * 1000 * (H0 / 3.086e19) / BETA  # m/s²
a0_kms2_kpc = a0_SI * 3.086e13  # (km/s)²/kpc

print("=" * 80)
print("  ATTEMPTING TO DERIVE ALL CBT PARAMETERS")
print("=" * 80)
print()
print(f"Fundamental: β = 2e = {BETA:.4f}")
print(f"Derived:     a₀ = cH₀/(2e) = {a0_SI:.3e} m/s²")
print()

# =============================================================================
# CANDIDATE DERIVATIONS FOR α₀
# =============================================================================

print("=" * 80)
print("  CANDIDATE DERIVATIONS FOR α₀")
print("=" * 80)
print()

alpha_candidates = {
    "1/β = 1/(2e)": 1/BETA,
    "1/√β = 1/√(2e)": 1/np.sqrt(BETA),
    "√(2/β) = √(1/e)": np.sqrt(1/e),
    "1/e": 1/e,
    "2/β = 1/e": 2/BETA,
    "1/π": 1/pi,
    "e/π²": e/pi**2,
    "√2/e (virial×entropy)": np.sqrt(2)/e,
    "1/2 (empirical)": 0.5,
}

for name, value in alpha_candidates.items():
    print(f"  {name:<30s} = {value:.4f}")

# =============================================================================
# CANDIDATE DERIVATIONS FOR s (log slope)
# =============================================================================

print()
print("=" * 80)
print("  CANDIDATE DERIVATIONS FOR s (log slope)")
print("=" * 80)
print()

# Physical reasoning for s:
# - Information scales as log(R) (holographic principle)
# - The slope should be related to dimensionality or fundamental ratios

s_candidates = {
    "1/e": 1/e,
    "1/π": 1/pi,
    "1/(2e)": 1/(2*e),
    "1/β = 1/(2e)": 1/BETA,
    "ln(2)/e": np.log(2)/e,
    "1/3 (spatial dimensions)": 1/3,
    "3/10 (D/R₀)": 3/10,
    "e/(π²)": e/pi**2,
    "1/4": 0.25,
    "0.30 (empirical)": 0.30,
}

for name, value in s_candidates.items():
    print(f"  {name:<30s} = {value:.4f}")

# =============================================================================
# CANDIDATE DERIVATIONS FOR r_th
# =============================================================================

print()
print("=" * 80)
print("  CANDIDATE DERIVATIONS FOR r_th")
print("=" * 80)
print()

# Current: r_th = 0.10 R + 2.0 kpc
# Physical reasoning:
# - Inner region where binding develops
# - Could be related to MOND transition radius
# - The 2.0 kpc could be a characteristic core size

print("Current formula: r_th = 0.10 R + 2.0 kpc")
print()
print("Possible interpretations:")
print("  • 0.10 = 1/10 (decimal scaling)")
print("  • 0.10 ≈ 1/(πe) = ", 1/(pi*e))
print("  • 0.10 ≈ 1/(2π) = ", 1/(2*pi))
print("  • 0.10 ≈ e/π² - 0.18 = ", e/pi**2)
print()
print("  • 2.0 kpc ≈ characteristic core radius")
print("  • 2.0 ≈ e/e = 1... × 2")
print("  • 2.0 ≈ π/e × something?", pi/e)
print()

# =============================================================================
# LOAD DATA AND TEST COMBINATIONS
# =============================================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "SPARC_data")

def load_sparc_galaxy(filepath):
    data = {'radius': [], 'velocity': [], 'error': [], 'v_bar': []}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 6:
                r = float(parts[0])
                v_obs = float(parts[1])
                v_err = float(parts[2])
                v_gas = float(parts[3])
                v_disk = float(parts[4])
                v_bul = float(parts[5])
                v_bar = np.sqrt(v_gas**2 + v_disk**2 + v_bul**2)
                if r > 0 and v_obs > 0 and v_err > 0:
                    data['radius'].append(r)
                    data['velocity'].append(v_obs)
                    data['error'].append(max(v_err, 1.0))
                    data['v_bar'].append(v_bar)
    if len(data['radius']) < 5:
        return None
    r_arr = np.array(data['radius'])
    v_arr = np.array(data['velocity'])
    return {
        'radius': r_arr,
        'velocity': v_arr,
        'error': np.array(data['error']),
        'v_bar': np.array(data['v_bar']),
        'v_max': max(v_arr),
        'size': max(r_arr),
    }

files = glob.glob(os.path.join(DATA_DIR, "*_rotmod.dat"))
galaxies = []
for f in files:
    data = load_sparc_galaxy(f)
    if data:
        galaxies.append(data)

print()
print(f"Testing on {len(galaxies)} galaxies...")
print()

def fit_galaxy(data, alpha_0, s, r_th_coef, r_th_offset):
    """Fit with given parameters."""
    r = data['radius']
    v = data['velocity']
    err = data['error']
    v_bar = data['v_bar']
    size = data['size']
    v_max = data['v_max']
    
    geometric = 1 + s * np.log10(max(size, 1) / 10)
    alpha = min(alpha_0 * geometric, 1.0)
    r_th = r_th_coef * size + r_th_offset
    
    try:
        def cbt_func(r_in, scale):
            v_bar_s = v_bar * scale
            v0_max = alpha * v_max
            v0 = v0_max * np.minimum(r_in / r_th, 1.0)
            return np.sqrt(v_bar_s**2 + v0**2)
        
        popt, _ = curve_fit(cbt_func, r, v, p0=[1.0], sigma=err,
                           bounds=([0.1], [3.0]), maxfev=5000)
        chi_cbt = np.sum(((v - cbt_func(r, *popt)) / err)**2) / (len(r) - 1)
    except:
        chi_cbt = float('inf')
    
    try:
        def newton_func(r_in, scale):
            return v_bar * scale
        popt_n, _ = curve_fit(newton_func, r, v, p0=[1.0], sigma=err,
                             bounds=([0.1], [3.0]), maxfev=5000)
        chi_n = np.sum(((v - newton_func(r, *popt_n)) / err)**2) / (len(r) - 1)
    except:
        chi_n = float('inf')
    
    return chi_cbt, chi_n

def evaluate(alpha_0, s, r_th_coef, r_th_offset):
    wins = 0
    total = 0
    chis = []
    for data in galaxies:
        chi_cbt, chi_n = fit_galaxy(data, alpha_0, s, r_th_coef, r_th_offset)
        if chi_cbt < 100 and chi_n < 100:
            if chi_cbt < chi_n:
                wins += 1
            total += 1
            chis.append(chi_cbt)
    return wins, total, np.mean(chis) if chis else float('inf')

# =============================================================================
# TEST DERIVED PARAMETER COMBINATIONS
# =============================================================================

print("=" * 80)
print("  TESTING DERIVED PARAMETER COMBINATIONS")
print("=" * 80)
print()

# Current empirical baseline
wins, total, chi = evaluate(0.50, 0.30, 0.10, 2.0)
print(f"EMPIRICAL BASELINE: α₀=0.50, s=0.30, r_th=0.10R+2")
print(f"  Wins: {wins}/{total} ({100*wins/total:.1f}%), χ² = {chi:.2f}")
print()

# Test derived combinations
derived_tests = [
    # (α₀, s, r_th_coef, r_th_offset, description)
    (1/np.sqrt(BETA), 1/e, 0.10, 2.0, "α₀=1/√β, s=1/e, r_th=empirical"),
    (1/np.sqrt(BETA), 1/3, 0.10, 2.0, "α₀=1/√β, s=1/3, r_th=empirical"),
    (1/np.sqrt(BETA), 0.30, 0.10, 2.0, "α₀=1/√β, s=0.30, r_th=empirical"),
    (np.sqrt(2)/e, 1/e, 0.10, 2.0, "α₀=√2/e, s=1/e, r_th=empirical"),
    (np.sqrt(2)/e, 1/3, 0.10, 2.0, "α₀=√2/e, s=1/3, r_th=empirical"),
    (1/e, 1/e, 0.10, 2.0, "α₀=1/e, s=1/e, r_th=empirical"),
    (1/e, 1/3, 0.10, 2.0, "α₀=1/e, s=1/3, r_th=empirical"),
    (1/e, 1/pi, 0.10, 2.0, "α₀=1/e, s=1/π, r_th=empirical"),
    # Try different r_th
    (1/np.sqrt(BETA), 1/e, 1/(2*pi), 2.0, "α₀=1/√β, s=1/e, r_th=(1/2π)R+2"),
    (1/np.sqrt(BETA), 1/e, 1/e, 2.0, "α₀=1/√β, s=1/e, r_th=(1/e)R+2"),
    (1/np.sqrt(BETA), 1/e, 0.10, e, "α₀=1/√β, s=1/e, r_th=0.1R+e"),
    (np.sqrt(2)/e, 1/e, 1/(2*pi), pi/e, "α₀=√2/e, s=1/e, r_th=(1/2π)R+π/e"),
]

print(f"{'Description':<50s} {'Wins':>8s} {'%':>6s} {'χ²':>8s}")
print("-" * 75)

for alpha_0, s, r_th_c, r_th_o, desc in derived_tests:
    wins, total, chi = evaluate(alpha_0, s, r_th_c, r_th_o)
    pct = 100*wins/total if total > 0 else 0
    print(f"{desc:<50s} {wins:>8d} {pct:>5.1f}% {chi:>8.2f}")

# =============================================================================
# FIND BEST FULLY-DERIVED FORMULA
# =============================================================================

print()
print("=" * 80)
print("  SEARCHING FOR BEST FULLY-DERIVED FORMULA")
print("=" * 80)
print()

# Candidate α₀ values derived from e
alpha_opts = [1/e, 1/np.sqrt(BETA), np.sqrt(2)/e, 2/BETA, 1/np.sqrt(e)]

# Candidate s values derived from fundamentals
s_opts = [1/e, 1/3, 1/pi, np.log(2)/e, e/pi**2, 1/4]

# Candidate r_th coefficients
rth_c_opts = [1/10, 1/(2*pi), 1/e, 1/pi]

# Candidate r_th offsets (in kpc, related to characteristic scales)
rth_o_opts = [2.0, e, pi/e, np.sqrt(2)*e, pi]

best_result = None
best_wins = 0
best_params = None

print("Testing all combinations of derived values...")
total_combos = len(alpha_opts) * len(s_opts) * len(rth_c_opts) * len(rth_o_opts)
print(f"({total_combos} combinations)")
print()

for a in alpha_opts:
    for s in s_opts:
        for rc in rth_c_opts:
            for ro in rth_o_opts:
                wins, total, chi = evaluate(a, s, rc, ro)
                if wins > best_wins:
                    best_wins = wins
                    best_params = (a, s, rc, ro)
                    best_result = (wins, total, chi)

a, s, rc, ro = best_params
wins, total, chi = best_result

print(f"BEST FULLY-DERIVED FORMULA:")
print(f"  α₀ = {a:.4f}")
print(f"  s = {s:.4f}")
print(f"  r_th = {rc:.4f} R + {ro:.4f}")
print()
print(f"  Wins: {wins}/{total} ({100*wins/total:.1f}%)")
print(f"  Mean χ²: {chi:.2f}")
print()

# Identify what these values are
def identify_value(val, name):
    candidates = {
        "1/e": 1/e,
        "1/√(2e)": 1/np.sqrt(BETA),
        "√2/e": np.sqrt(2)/e,
        "1/√e": 1/np.sqrt(e),
        "2/(2e)": 2/BETA,
        "1/3": 1/3,
        "1/π": 1/pi,
        "ln(2)/e": np.log(2)/e,
        "e/π²": e/pi**2,
        "1/4": 0.25,
        "1/10": 0.1,
        "1/(2π)": 1/(2*pi),
        "e": e,
        "π/e": pi/e,
        "√2·e": np.sqrt(2)*e,
        "π": pi,
        "2": 2.0,
    }
    for cname, cval in candidates.items():
        if abs(val - cval) < 0.001:
            print(f"  {name} = {cname}")
            return
    print(f"  {name} = {val:.4f} (no simple form found)")

print("Identifying best parameters:")
identify_value(a, "α₀")
identify_value(s, "s")
identify_value(rc, "r_th coefficient")
identify_value(ro, "r_th offset")

print()
print("=" * 80)
print("  CONCLUSION")
print("=" * 80)
