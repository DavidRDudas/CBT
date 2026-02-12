#!/usr/bin/env python3
"""
CBT CMB Comprehensive Test: TT + TE + EE vs Planck 2018
========================================================

Tests CBT binding field prediction against ALL three main Planck 
CMB power spectra: TT, TE, and EE.

CBT prediction: omega_scf = omega_b × 2e = 0.1216 (fixed)
Remaining 5 params: h=0.666, A_s=2.1e-9, n_s=0.960, tau=0.052

Author: D. Dudas, 2026
"""

import numpy as np
import subprocess
import os
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CLASS_DIR = os.path.join(SCRIPT_DIR, "class")
CLASS_BIN = os.path.join(CLASS_DIR, "class")
OUTPUT_DIR = os.path.join(CLASS_DIR, "output")
DATA_DIR = SCRIPT_DIR

omega_b = 0.02237
omega_cdm_cbt = omega_b * 2 * np.e  # 0.1216

# ============================================================
# CLASS column mapping:
# col 0: l, col 1: TT, col 2: EE, col 3: TE, col 4: BB, col 5: phiphi
# All are dimensionless l(l+1)C_l/(2pi)
# ============================================================

models = {
    'lcdm': {
        'label': r'$\Lambda$CDM (6 params)',
        'h': 0.6736, 'omega_cdm': 0.1200, 'A_s': 2.1e-9, 
        'n_s': 0.9649, 'tau': 0.0544,
        'color': '#2196F3', 'ls': '-', 'lw': 2.0,
    },
    'cbt': {
        'label': r'CBT (5 params, $\omega_{cdm}$ predicted)',
        'h': 0.666, 'omega_cdm': omega_cdm_cbt, 'A_s': 2.1e-9,
        'n_s': 0.960, 'tau': 0.052,
        'color': '#FF5722', 'ls': '--', 'lw': 1.8,
    },
}

print("=" * 65)
print("  CBT CMB COMPREHENSIVE TEST: TT + TE + EE")  
print("=" * 65)

# ============================================================
# Run CLASS
# ============================================================
T_cmb_uK = 2.7255e6
spectra = {}

for key, m in models.items():
    label = f"cmb4_{key}"
    ini = f"""
output = tCl,pCl,lCl
lensing = yes
h = {m['h']}
T_cmb = 2.7255
omega_b = {omega_b}
omega_cdm = {m['omega_cdm']}
N_ur = 2.0328
N_ncdm = 1
m_ncdm = 0.06
A_s = {m['A_s']}
n_s = {m['n_s']}
tau_reio = {m['tau']}
k_pivot = 0.05
l_max_scalars = 2500
format = class
headers = yes
root = output/{label}_
background_verbose = 0
thermodynamics_verbose = 0
perturbations_verbose = 0
transfer_verbose = 0
primordial_verbose = 0
harmonic_verbose = 0
lensing_verbose = 0
output_verbose = 0
"""
    ini_path = os.path.join(CLASS_DIR, f"_temp_{label}.ini")
    with open(ini_path, 'w') as f:
        f.write(ini)
    
    print(f"\n  Running CLASS: {key}... ", end="", flush=True)
    result = subprocess.run([CLASS_BIN, ini_path], capture_output=True, cwd=CLASS_DIR)
    
    if result.returncode != 0:
        print("FAILED")
        continue
    print("✓")
    
    data = np.loadtxt(os.path.join(OUTPUT_DIR, f"{label}_00_cl_lensed.dat"))
    ell = data[:, 0]
    spectra[key] = {
        'ell': ell,
        'TT': data[:, 1] * T_cmb_uK**2,             # μK²
        'EE': data[:, 2] * T_cmb_uK**2,             # μK²
        'TE': data[:, 3] * T_cmb_uK**2,             # μK²
    }

# ============================================================
# Load Planck data for TT, TE, EE
# ============================================================
planck = {}

# TT
tt_data = np.loadtxt(os.path.join(DATA_DIR, "COM_PowerSpect_CMB-TT-binned_R3.01.txt"))
planck['TT'] = {
    'ell': tt_data[:, 0], 
    'Dl': tt_data[:, 1], 
    'err': (tt_data[:, 2] + tt_data[:, 3]) / 2,
    'bestfit': tt_data[:, 4]
}

# TE
te_data = np.loadtxt(os.path.join(DATA_DIR, "COM_PowerSpect_CMB-TE-binned_R3.02.txt"))
planck['TE'] = {
    'ell': te_data[:, 0],
    'Dl': te_data[:, 1],
    'err': (te_data[:, 2] + te_data[:, 3]) / 2,
    'bestfit': te_data[:, 4]
}

# EE  
ee_data = np.loadtxt(os.path.join(DATA_DIR, "COM_PowerSpect_CMB-EE-binned_R3.02.txt"))
planck['EE'] = {
    'ell': ee_data[:, 0],
    'Dl': ee_data[:, 1],
    'err': (ee_data[:, 2] + ee_data[:, 3]) / 2,
    'bestfit': ee_data[:, 4]
}

# ============================================================
# Compute χ² for each spectrum
# ============================================================
print()
print("=" * 65)
print("  RESULTS: χ²/dof by spectrum")
print("=" * 65)
print()
print(f"  {'Spectrum':<10s}  {'ΛCDM':>12s}  {'CBT':>12s}  {'Planck bins':>12s}")
print(f"  {'-'*50}")

total_chi2 = {'lcdm': 0, 'cbt': 0}
total_dof = 0

for spec_name in ['TT', 'TE', 'EE']:
    p = planck[spec_name]
    n_bins = len(p['ell'])
    dof = n_bins - 1
    total_dof += dof
    
    chi2_vals = {}
    for key in ['lcdm', 'cbt']:
        s = spectra[key]
        f_model = interp1d(s['ell'], s[spec_name], bounds_error=False, fill_value=0)
        dl_pred = f_model(p['ell'])
        chi2 = float(np.sum(((p['Dl'] - dl_pred) / p['err'])**2))
        chi2_vals[key] = chi2
        total_chi2[key] += chi2
    
    print(f"  {spec_name:<10s}  {chi2_vals['lcdm']/dof:>10.3f}    {chi2_vals['cbt']/dof:>10.3f}    {n_bins:>8d}")

print(f"  {'-'*50}")
print(f"  {'COMBINED':<10s}  {total_chi2['lcdm']/total_dof:>10.3f}    {total_chi2['cbt']/total_dof:>10.3f}    {total_dof+3:>8d}")

print()
print(f"  ΛCDM combined χ²/dof = {total_chi2['lcdm']/total_dof:.4f}")
print(f"  CBT  combined χ²/dof = {total_chi2['cbt']/total_dof:.4f}")
print()

if total_chi2['cbt']/total_dof < 1.5:
    print(f"  ✅ CBT PASSES THE COMPREHENSIVE CMB TEST")
elif total_chi2['cbt']/total_dof < 2.0:
    print(f"  ⚠️  CBT IS MARGINALLY CONSISTENT")
else:
    print(f"  ❌ CBT FAILS THE COMPREHENSIVE CMB TEST")

# ============================================================
# Publication-quality 3-panel plot
# ============================================================
fig, axes = plt.subplots(3, 2, figsize=(16, 14),
                          gridspec_kw={'width_ratios': [3, 1], 'hspace': 0.25, 'wspace': 0.05})
plt.rcParams.update({'font.family': 'serif', 'font.size': 12})

spec_configs = {
    'TT': {'ylabel': r'$\mathcal{D}_\ell^{TT}$ [$\mu$K$^2$]', 'row': 0, 'ylim_main': (-200, 6500)},
    'TE': {'ylabel': r'$\mathcal{D}_\ell^{TE}$ [$\mu$K$^2$]', 'row': 1, 'ylim_main': (-160, 120)},
    'EE': {'ylabel': r'$\mathcal{D}_\ell^{EE}$ [$\mu$K$^2$]', 'row': 2, 'ylim_main': (-5, 45)},
}

for spec_name, cfg in spec_configs.items():
    ax_main = axes[cfg['row'], 0]
    ax_res = axes[cfg['row'], 1]
    p = planck[spec_name]
    
    # Main panel
    ax_main.errorbar(p['ell'], p['Dl'], yerr=p['err'], fmt='o', color='black',
                     markersize=2, alpha=0.4, zorder=5, capsize=0, elinewidth=0.5,
                     label='Planck 2018')
    
    for key in ['lcdm', 'cbt']:
        m = models[key]
        s = spectra[key]
        f_model = interp1d(s['ell'], s[spec_name], bounds_error=False, fill_value=0)
        chi2 = float(np.sum(((p['Dl'] - f_model(p['ell'])) / p['err'])**2))
        dof = len(p['ell']) - 1
        ax_main.plot(s['ell'], s[spec_name], ls=m['ls'], color=m['color'], 
                     lw=m['lw'], alpha=0.9,
                     label=f"{m['label']} (χ²/dof={chi2/dof:.2f})")
    
    ax_main.set_ylabel(cfg['ylabel'], fontsize=13)
    ax_main.set_xlim(2, 2500)
    ax_main.set_ylim(cfg['ylim_main'])
    ax_main.grid(True, alpha=0.15)
    ax_main.legend(fontsize=9, loc='upper right' if spec_name == 'TT' else 'best')
    
    if spec_name != 'EE':
        ax_main.set_xticklabels([])
    else:
        ax_main.set_xlabel(r'Multipole $\ell$', fontsize=14)
    
    # Residual: CBT vs LCDM
    s_l = spectra['lcdm']
    s_c = spectra['cbt']
    f_c = interp1d(s_c['ell'], s_c[spec_name], bounds_error=False, fill_value=0)
    
    if spec_name == 'TT':
        res = (f_c(s_l['ell']) - s_l[spec_name]) / s_l[spec_name] * 100
        ax_res.set_ylabel('Δ/D [%]', fontsize=10)
        ax_res.set_ylim(-5, 5)
    else:
        # For TE/EE, use absolute residual
        res = f_c(s_l['ell']) - s_l[spec_name]
        ax_res.set_ylabel(r'Δ [$\mu$K$^2$]', fontsize=10)
    
    ax_res.plot(s_l['ell'], res, '-', color='#FF5722', lw=1.0, alpha=0.7)
    ax_res.axhline(0, color='gray', ls='--', alpha=0.5)
    ax_res.set_xlim(2, 2500)
    ax_res.grid(True, alpha=0.15)
    ax_res.set_title(f'{spec_name} residual', fontsize=10)
    
    if spec_name != 'EE':
        ax_res.set_xticklabels([])
    else:
        ax_res.set_xlabel(r'$\ell$', fontsize=12)

fig.suptitle('CMB Power Spectra: CBT Binding Field vs ΛCDM vs Planck 2018', 
             fontsize=16, fontweight='bold', y=0.98)

# Text box
combined_lcdm = total_chi2['lcdm']/total_dof
combined_cbt = total_chi2['cbt']/total_dof
textstr = (f'Combined TT+TE+EE:\n'
           f'  ΛCDM: χ²/dof = {combined_lcdm:.3f} (6 params)\n'
           f'  CBT:  χ²/dof = {combined_cbt:.3f} (5 params)\n'
           f'  ω_cdm = ω_b×2e = {omega_cdm_cbt:.4f}')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
axes[0, 0].text(0.02, 0.97, textstr, transform=axes[0, 0].transAxes, fontsize=10,
                verticalalignment='top', bbox=props, family='monospace')

plt.tight_layout(rect=[0, 0, 1, 0.96])

out_path = os.path.join(DATA_DIR, "cbt_cmb_all_spectra.png")
plt.savefig(out_path, dpi=200, bbox_inches='tight')
print(f"\n✓ Plot saved: {out_path}")
plt.close()

print()
print("=" * 65)
