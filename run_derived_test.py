#!/usr/bin/env python3
"""
CBT CANONICAL TEST - DERIVED FORMULA VERSION
=============================================

This version uses the fully-derived parameters:
- α₀ = 1/e ≈ 0.368
- s = 1/e ≈ 0.368
- r_th = R/(2π) + √2·e ≈ 0.159R + 3.85 kpc

All parameters derived from e and π.

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
# DERIVED CONSTANTS
# =============================================================================

E = np.e
PI = np.pi
BETA = 2 * E

# All derived from e and π:
ALPHA_0 = 1 / E        # ≈ 0.368
S = 1 / E              # ≈ 0.368
R0 = 10                # kpc (reference scale)
R_TH_COEF = 1 / (2*PI) # ≈ 0.159
R_TH_OFFSET = np.sqrt(2) * E  # ≈ 3.85 kpc

# MOND acceleration scale: a₀ = cH₀/(2e)
# In units of (km/s)²/kpc for direct use with SPARC data
A0 = 2.998e5 * (70 / 3.086e19) / BETA * 3.086e16  # ≈ 3480 (km/s)²/kpc

print(f"Derived parameters from e = {E:.4f} and π = {PI:.4f}:")
print(f"  α₀ = 1/e = {ALPHA_0:.4f}")
print(f"  s = 1/e = {S:.4f}")
print(f"  r_th = R/{2*PI:.4f} + {R_TH_OFFSET:.4f} kpc")
print()

# =============================================================================
# CONFIGURATION
# =============================================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "SPARC_data")

# =============================================================================
# DATA LOADING  
# =============================================================================

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
    
    return {
        'radius': np.array(data['radius']),
        'velocity': np.array(data['velocity']),
        'error': np.array(data['error']),
        'v_bar': np.array(data['v_bar']),
        'v_max': max(data['velocity']),
        'size': max(data['radius'])
    }

# =============================================================================
# DERIVED MODEL
# =============================================================================

def alpha_derived(R_kpc):
    """
    Binding strength - FULLY DERIVED from thermodynamics.
    α(R) = (1/e) × (1 + (1/e) × log₁₀(R/10))
    """
    return min(ALPHA_0 * (1 + S * np.log10(max(R_kpc, 1) / R0)), 1.0)


def cbt_model_derived(r, v_bar, v_max, size, scale):
    """
    CBT model with derived parameters.
    """
    v_bar_scaled = v_bar * scale
    
    alpha = alpha_derived(size)
    r_th = R_TH_COEF * size + R_TH_OFFSET
    
    v0_max = alpha * v_max
    v0 = v0_max * np.minimum(r / r_th, 1.0)
    
    return np.sqrt(v_bar_scaled**2 + v0**2)


def newtonian_model(r, v_bar, scale):
    return v_bar * scale


def mond_model(r, v_bar, scale):
    """
    Standard MOND with μ(x) = x/√(1+x²).
    g_obs = g_N / μ(g_N/a₀), v = √(g_obs × r)
    """
    v_bar_scaled = v_bar * scale
    g_N = np.maximum(v_bar_scaled**2 / (r + 0.01), 1e-20)
    x = g_N / A0
    mu = x / np.sqrt(1 + x**2)
    g_obs = g_N / mu
    return np.sqrt(np.maximum(g_obs * r, 0))


def fit_galaxy(data):
    """Fit Newton, CBT, and MOND to a galaxy."""
    r = data['radius']
    v = data['velocity']
    err = data['error']
    v_bar = data['v_bar']
    size = data['size']
    v_max = data['v_max']
    
    # CBT fit (1 free parameter: M/L scale)
    try:
        def cbt_func(r_in, scale):
            return cbt_model_derived(r_in, v_bar, v_max, size, scale)
        
        popt_c, _ = curve_fit(cbt_func, r, v, p0=[1.0], sigma=err,
                             bounds=([0.1], [3.0]), maxfev=5000)
        v_pred_c = cbt_func(r, *popt_c)
        chi_c = np.sum(((v - v_pred_c) / err)**2) / (len(r) - 1)
    except:
        chi_c = float('inf')
    
    # Newton fit (1 free parameter: M/L scale)
    try:
        def newton_func(r_in, scale):
            return newtonian_model(r_in, v_bar, scale)
        
        popt_n, _ = curve_fit(newton_func, r, v, p0=[1.0], sigma=err,
                             bounds=([0.1], [3.0]), maxfev=5000)
        v_pred_n = newton_func(r, *popt_n)
        chi_n = np.sum(((v - v_pred_n) / err)**2) / (len(r) - 1)
    except:
        chi_n = float('inf')
    
    # MOND fit (1 free parameter: M/L scale)
    try:
        def mond_func(r_in, scale):
            return mond_model(r_in, v_bar, scale)
        
        popt_m, _ = curve_fit(mond_func, r, v, p0=[1.0], sigma=err,
                             bounds=([0.1], [3.0]), maxfev=5000)
        v_pred_m = mond_func(r, *popt_m)
        chi_m = np.sum(((v - v_pred_m) / err)**2) / (len(r) - 1)
    except:
        chi_m = float('inf')
    
    return chi_c, chi_n, chi_m, alpha_derived(size)

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("  COMPLEXITY BINDING THEORY - DERIVED FORMULA VALIDATION")
    print("=" * 70)
    print()
    print("  All parameters derived from e and π:")
    print(f"    α₀ = 1/e = {ALPHA_0:.4f}")
    print(f"    s = 1/e = {S:.4f}")
    print(f"    r_th = R/(2π) + √2·e = {R_TH_COEF:.4f}R + {R_TH_OFFSET:.2f} kpc")
    print()
    
    # Load galaxies
    if not os.path.exists(DATA_DIR):
        print(f"ERROR: Data directory not found: {DATA_DIR}")
        return
    
    files = glob.glob(os.path.join(DATA_DIR, "*_rotmod.dat"))
    print(f"  Found {len(files)} galaxy files")
    print()
    
    # Evaluate
    results = []
    cbt_wins_N = 0
    newton_wins = 0
    cbt_wins_M = 0
    mond_wins = 0
    
    for filepath in sorted(files):
        name = os.path.basename(filepath).replace("_rotmod.dat", "")
        data = load_sparc_galaxy(filepath)
        
        if data is None:
            continue
        
        chi_c, chi_n, chi_m, alpha = fit_galaxy(data)
        
        if chi_c == float('inf') or chi_n == float('inf'):
            continue
        if chi_c > 100 or chi_n > 100:
            continue
        
        if chi_c < chi_n:
            cbt_wins_N += 1
            winner = "CBT"
        else:
            newton_wins += 1
            winner = "Newton"
        
        if chi_m < float('inf') and chi_m < 100:
            if chi_c < chi_m:
                cbt_wins_M += 1
            else:
                mond_wins += 1
        
        results.append({
            'name': name,
            'chi_cbt': chi_c,
            'chi_newton': chi_n,
            'chi_mond': chi_m,
            'winner': winner,
            'alpha': alpha,
            'v_max': data['v_max'],
            'size': data['size']
        })
    
    total = len(results)
    
    print("=" * 70)
    print("  RESULTS")
    print("=" * 70)
    print()
    mond_valid = cbt_wins_M + mond_wins
    
    print(f"  Valid galaxies: {total}")
    print()
    print(f"  --- CBT vs Newton ---")
    print(f"  CBT wins:    {cbt_wins_N} ({100*cbt_wins_N/total:.1f}%)")
    print(f"  Newton wins: {newton_wins} ({100*newton_wins/total:.1f}%)")
    print()
    if mond_valid > 0:
        print(f"  --- CBT vs MOND ---")
        print(f"  CBT wins:    {cbt_wins_M} ({100*cbt_wins_M/mond_valid:.1f}%)")
        print(f"  MOND wins:   {mond_wins} ({100*mond_wins/mond_valid:.1f}%)")
        print()
    
    avg_cbt = np.mean([r['chi_cbt'] for r in results])
    avg_newton = np.mean([r['chi_newton'] for r in results])
    mond_v = [r for r in results if r['chi_mond'] < 100]
    avg_mond = np.mean([r['chi_mond'] for r in mond_v]) if mond_v else float('nan')
    
    print(f"  Mean χ² (Newton): {avg_newton:.2f}")
    print(f"  Mean χ² (CBT):    {avg_cbt:.2f}")
    print(f"  Mean χ² (MOND):   {avg_mond:.2f}")
    print(f"  CBT improvement over Newton: {avg_n/avg_c:.1f}x" if 'avg_n' in dir() else "")
    print()
    print("=" * 70)
    
    # Save results
    output_file = os.path.join(SCRIPT_DIR, "results_derived_formula.csv")
    with open(output_file, 'w') as f:
        f.write("galaxy,chi_cbt,chi_newton,chi_mond,winner,alpha,v_max,size\n")
        for r in results:
            f.write(f"{r['name']},{r['chi_cbt']:.4f},{r['chi_newton']:.4f},"
                    f"{r.get('chi_mond', 'inf'):.4f},"
                    f"{r['winner']},{r['alpha']:.3f},{r['v_max']:.1f},{r['size']:.1f}\n")
    print(f"  Results saved to: {output_file}")
    print("=" * 70)


if __name__ == "__main__":
    main()
