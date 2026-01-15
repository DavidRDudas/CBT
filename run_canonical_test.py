"""
This is the official, reproducible implementation of Complexity Binding Theory
testing on the SPARC database.

Usage:
    python run_canonical_test.py [--use-btf-vmax]

Author: David R. Dudas
"""

import numpy as np
from scipy.optimize import curve_fit
import os
import glob
import argparse
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

# Use relative path - expects SPARC_data folder in same directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "SPARC_data")

# Gravitational constant in useful units
G = 4.3e-6  # kpc * (km/s)^2 / M_sun

# =============================================================================
# DATA LOADING
# =============================================================================

def load_sparc_galaxy(filepath):
    """Load a SPARC rotation curve file."""
    data = {'radius': [], 'velocity': [], 'error': [], 'v_bar': []}
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 6:
                r = float(parts[0])       # kpc
                v_obs = float(parts[1])   # km/s
                v_err = float(parts[2])   # km/s
                v_gas = float(parts[3])   # km/s
                v_disk = float(parts[4])  # km/s
                v_bul = float(parts[5])   # km/s
                
                # Baryonic velocity
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
# MODEL DEFINITIONS
# =============================================================================

def alpha_universal(R_kpc):
    """
    Universal binding strength formula.
    This is FIXED, not fitted per galaxy.
    """
    return min(0.50 * (1 + 0.3 * np.log10(max(R_kpc, 1) / 10)), 1.0)


def newtonian_model(r, v_bar, scale):
    """Newtonian prediction: just scaled baryonic matter."""
    return v_bar * scale


def cbt_model(r, v_bar, v_max, size, scale):
    """
    Complexity Binding Theory model.
    
    Parameters:
        r: radius array
        v_bar: baryonic velocity array  
        v_max: maximum observed velocity (or BTF-predicted)
        size: galaxy size in kpc
        scale: M/L scaling factor (the ONLY fitted parameter)
    
    The α(R) is computed from the universal formula, NOT fitted.
    """
    # Apply M/L scaling to baryonic component
    v_bar_scaled = v_bar * scale
    
    # Universal α from formula (NOT a free parameter)
    alpha = alpha_universal(size)
    
    # Binding velocity
    v0_max = alpha * v_max
    r_threshold = 0.10 * size + 2.0
    v0 = v0_max * np.minimum(r / r_threshold, 1.0)
    
    # Total velocity
    return np.sqrt(v_bar_scaled**2 + v0**2)

# =============================================================================
# FITTING FUNCTIONS
# =============================================================================

def fit_galaxy(data, use_btf_vmax=False):
    """
    Fit both Newton and CBT models to a galaxy.
    
    Both models have exactly 1 free parameter: the M/L scale.
    CBT uses the universal α(R) formula with no additional freedom.
    """
    r = data['radius']
    v = data['velocity']
    err = data['error']
    v_bar = data['v_bar']
    size = data['size']
    
    # V_max: use observed or BTF-predicted
    if use_btf_vmax:
        # Baryonic Tully-Fisher: M_bar = A * V_flat^4
        # Approximate V_flat from baryonic velocity at large r
        v_bar_max = np.max(v_bar)
        v_max = v_bar_max * 1.5  # Rough BTF-based estimate
    else:
        v_max = data['v_max']
    
    # Newtonian fit (1 parameter: scale)
    try:
        def newton_func(r_in, scale):
            return newtonian_model(r_in, v_bar, scale)
        
        popt_n, _ = curve_fit(newton_func, r, v, p0=[1.0], sigma=err,
                             bounds=([0.1], [3.0]), maxfev=5000)
        v_pred_n = newton_func(r, *popt_n)
        chi_n = np.sum(((v - v_pred_n) / err)**2) / (len(r) - 1)
        scale_n = popt_n[0]
    except:
        chi_n = float('inf')
        scale_n = 1.0
    
    # CBT fit (1 parameter: scale; α is universal, not fitted)
    try:
        def cbt_func(r_in, scale):
            return cbt_model(r_in, v_bar, v_max, size, scale)
        
        popt_c, _ = curve_fit(cbt_func, r, v, p0=[1.0], sigma=err,
                             bounds=([0.1], [3.0]), maxfev=5000)
        v_pred_c = cbt_func(r, *popt_c)
        chi_c = np.sum(((v - v_pred_c) / err)**2) / (len(r) - 1)
        scale_c = popt_c[0]
    except:
        chi_c = float('inf')
        scale_c = 1.0
    
    return {
        'chi_newton': chi_n,
        'chi_cbt': chi_c,
        'scale_newton': scale_n,
        'scale_cbt': scale_c,
        'alpha_used': alpha_universal(size)
    }

# =============================================================================
# MAIN EVALUATION
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Canonical CBT SPARC Test')
    parser.add_argument('--use-btf-vmax', action='store_true',
                        help='Use BTF-predicted Vmax instead of observed')
    args = parser.parse_args()
    
    print("=" * 70)
    print("  COMPLEXITY BINDING THEORY - CANONICAL SPARC EVALUATION")
    print("=" * 70)
    print()
    print("  Parameter accounting:")
    print("    - Newton: 1 free parameter (M/L scale)")
    print("    - CBT:    1 free parameter (M/L scale) + universal α(R)")
    print()
    print(f"  V_max source: {'BTF-predicted' if args.use_btf_vmax else 'Observed'}")
    print()
    
    # Load galaxies
    if not os.path.exists(DATA_DIR):
        print(f"  ERROR: Data directory not found: {DATA_DIR}")
        print("  Please ensure SPARC_data folder is in the same directory.")
        return
    
    files = glob.glob(os.path.join(DATA_DIR, "*_rotmod.dat"))
    print(f"  Found {len(files)} galaxy files")
    print()
    
    # Run evaluation
    results = []
    newton_wins = 0
    cbt_wins = 0
    
    for filepath in sorted(files):
        name = os.path.basename(filepath).replace("_rotmod.dat", "")
        data = load_sparc_galaxy(filepath)
        
        if data is None:
            continue
        
        fit_result = fit_galaxy(data, use_btf_vmax=args.use_btf_vmax)
        
        chi_n = fit_result['chi_newton']
        chi_c = fit_result['chi_cbt']
        
        if chi_n == float('inf') and chi_c == float('inf'):
            continue
        
        if chi_c < chi_n:
            winner = "CBT"
            cbt_wins += 1
        else:
            winner = "Newton"
            newton_wins += 1
        
        results.append({
            'name': name,
            'chi_newton': chi_n,
            'chi_cbt': chi_c,
            'winner': winner,
            'alpha': fit_result['alpha_used'],
            'v_max': data['v_max'],
            'size': data['size']
        })
    
    # Summary
    total = len(results)
    print("=" * 70)
    print("  RESULTS")
    print("=" * 70)
    print()
    print(f"  Total galaxies: {total}")
    print()
    print(f"  CBT wins:    {cbt_wins} ({100*cbt_wins/total:.1f}%)")
    print(f"  Newton wins: {newton_wins} ({100*newton_wins/total:.1f}%)")
    print()
    
    # Average chi-squared
    valid = [r for r in results if r['chi_newton'] < 100 and r['chi_cbt'] < 100]
    avg_n = np.mean([r['chi_newton'] for r in valid])
    avg_c = np.mean([r['chi_cbt'] for r in valid])
    
    print(f"  Mean χ² (Newton): {avg_n:.2f}")
    print(f"  Mean χ² (CBT):    {avg_c:.2f}")
    print(f"  Improvement:      {avg_n/avg_c:.1f}x")
    print()
    print("=" * 70)
    
    # Save results to CSV
    output_file = os.path.join(SCRIPT_DIR, "results_canonical.csv")
    with open(output_file, 'w') as f:
        f.write("galaxy,chi_newton,chi_cbt,winner,alpha,v_max,size\n")
        for r in results:
            f.write(f"{r['name']},{r['chi_newton']:.4f},{r['chi_cbt']:.4f},"
                    f"{r['winner']},{r['alpha']:.3f},{r['v_max']:.1f},{r['size']:.1f}\n")
    print(f"  Results saved to: {output_file}")
    print("=" * 70)

if __name__ == "__main__":
    main()
