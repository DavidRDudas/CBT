#!/usr/bin/env python3
"""
CBT vs MOND vs Newton: Using EXACT Paper 1 methodology
========================================================

Both CBT and MOND: v² = v_bar_scaled² + v_extra²
  CBT:   v_extra = v₀ = α(R) × v_max × min(r/r_th, 1)
  MOND:  v_extra comes from μ(g/a₀) = x/√(1+x²)
         i.e. g_obs = g_N/μ(g_N/a₀), v_mond² = g_obs × r
  Newton: v = v_bar_scaled

All have 1 free parameter: M/L scale.

Trying with β=6 params (α₀=0.50, s=0.30, r_th=0.10R+2.0)
"""

import numpy as np
from scipy.optimize import curve_fit, minimize_scalar
import os, glob
import warnings
warnings.filterwarnings('ignore')

# β=6 params
ALPHA_0_B6 = 0.50; S_B6 = 0.30

# β=2e params  
ALPHA_0_D = 1/np.e; S_D = 1/np.e

R0 = 10
a0_b6 = 2.998e8 * (70/3.086e19) / 6.0 * 3.086e13
a0_b2e = 2.998e8 * (70/3.086e19) / (2*np.e) * 3.086e13

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "SPARC_data")
files = sorted(glob.glob(os.path.join(data_dir, "*_rotmod.dat")))

galaxies = []
for fpath in files:
    name = os.path.basename(fpath).replace("_rotmod.dat", "")
    rows = []
    with open(fpath, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.split()
            if len(p) >= 6:
                r,vo,ve,vg,vd,vb = [float(x) for x in p[:6]]
                vbar = np.sqrt(vg**2 + vd**2 + vb**2)
                if r > 0 and vo > 0 and ve > 0:
                    rows.append((r, vo, max(ve,1.0), vbar))
    if len(rows) < 5: continue
    d = np.array(rows)
    galaxies.append({
        'name': name, 'r': d[:,0], 'v': d[:,1], 'err': d[:,2], 'vbar': d[:,3],
        'vmax': max(d[:,1]), 'size': max(d[:,0])
    })

# ======================================================================
# Models
# ======================================================================
def cbt_v(r, vbar, vmax, size, scale, alpha0, s):
    vs = vbar * scale
    alpha = min(alpha0 * (1 + s * np.log10(max(size, 1) / R0)), 1.0)
    rth = 0.10 * size + 2.0
    v0 = alpha * vmax * np.minimum(r / rth, 1.0)
    return np.sqrt(vs**2 + v0**2)

def mond_v(r, vbar, scale, a0):
    """MOND with μ(x)=x/√(1+x²)"""
    vs = vbar * scale
    gN = np.maximum(vs**2 / (r + 0.01), 1e-20)
    x = gN / a0
    mu = x / np.sqrt(1 + x**2)
    gobs = gN / mu
    return np.sqrt(np.maximum(gobs * r, 0))

def newton_v(r, vbar, scale):
    return vbar * scale

# ======================================================================
# Fit and compare
# ======================================================================
def fit_model(r, v, err, model_func):
    """Fit a 1-parameter model (scale) using curve_fit."""
    try:
        popt, _ = curve_fit(model_func, r, v, p0=[1.0], sigma=err, 
                           bounds=([0.1],[3.0]), maxfev=5000)
        vpred = model_func(r, *popt)
        chi2 = np.sum(((v - vpred)/err)**2) / max(len(r)-1, 1)
        return chi2, popt[0]
    except:
        return float('inf'), 1.0

print("=" * 70)
print("  CBT vs MOND vs NEWTON — Paper 1 Methodology")
print("=" * 70)

for beta_label, a0, a0_val, s_val in [
    ('β=6', ALPHA_0_B6, a0_b6, S_B6),
    ('β=2e', ALPHA_0_D, a0_b2e, S_D),
]:
    results = []
    
    for gal in galaxies:
        r, v, err, vbar = gal['r'], gal['v'], gal['err'], gal['vbar']
        vmax, size = gal['vmax'], gal['size']
    
        chi_cbt, _ = fit_model(r, v, err, 
            lambda ri, s: cbt_v(ri, vbar, vmax, size, s, a0, s_val))
        chi_mond, _ = fit_model(r, v, err, 
            lambda ri, s: mond_v(ri, vbar, s, a0_val))
        chi_N, _ = fit_model(r, v, err, 
            lambda ri, s: newton_v(ri, vbar, s))
    
        results.append({'name': gal['name'], 
                       'chi_cbt': chi_cbt, 'chi_mond': chi_mond, 'chi_N': chi_N})

    n = len(results)
    cbt_vs_N = sum(1 for r in results if r['chi_cbt'] < r['chi_N'])
    cbt_vs_M = sum(1 for r in results if r['chi_cbt'] < r['chi_mond'])
    
    # Without filter
    print(f"\n  --- {beta_label} (all {n} galaxies, NO filter) ---")
    print(f"  CBT vs Newton: {cbt_vs_N}/{n} = {100*cbt_vs_N/n:.1f}%")
    print(f"  CBT vs MOND:   {cbt_vs_M}/{n} = {100*cbt_vs_M/n:.1f}%")
    
    # With chi2 < 100 filter
    filt = [r for r in results 
            if r['chi_cbt'] < 100 and r['chi_mond'] < 100 and r['chi_N'] < 100]
    nf = len(filt)
    cbt_vs_N_f = sum(1 for r in filt if r['chi_cbt'] < r['chi_N'])
    cbt_vs_M_f = sum(1 for r in filt if r['chi_cbt'] < r['chi_mond'])
    
    print(f"\n  --- {beta_label} ({nf} galaxies, χ²<100 filter) ---")
    print(f"  CBT vs Newton: {cbt_vs_N_f}/{nf} = {100*cbt_vs_N_f/nf:.1f}%")
    print(f"  CBT vs MOND:   {cbt_vs_M_f}/{nf} = {100*cbt_vs_M_f/nf:.1f}%")
    
    print(f"\n  Avg χ²/dof:")
    print(f"    CBT:  {np.mean([r['chi_cbt'] for r in filt]):.2f}")
    print(f"    MOND: {np.mean([r['chi_mond'] for r in filt]):.2f}")
    print(f"    Newton: {np.mean([r['chi_N'] for r in filt]):.2f}")
    
    # How many does MOND beat Newton?
    mond_vs_N = sum(1 for r in filt if r['chi_mond'] < r['chi_N'])
    print(f"\n  MOND vs Newton: {mond_vs_N}/{nf} = {100*mond_vs_N/nf:.1f}%")
    
    # Detailed wins/losses for CBT vs MOND
    big_wins = sorted([r for r in filt if r['chi_cbt'] < r['chi_mond']], 
                      key=lambda x: x['chi_mond']/max(x['chi_cbt'],0.01))
    print(f"\n  Top 10 CBT wins over MOND:")
    for r in big_wins[-10:][::-1]:
        print(f"    {r['name']:<18} CBT={r['chi_cbt']:.2f}  MOND={r['chi_mond']:.2f}  "
              f"N={r['chi_N']:.2f}  ratio={r['chi_mond']/max(r['chi_cbt'],0.01):.1f}x")
