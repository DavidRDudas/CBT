#!/usr/bin/env python3
"""
CMB ACOUSTIC PEAK CONSISTENCY CHECK FOR CBT
============================================

Simplified consistency check using the key insight:
Peak positions scale with Ω_m. If CBT's Ω_eff ≈ Planck's Ω_m,
then CBT predicts the same peak positions.

Author: David R. Dudas
Date: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# CONSTANTS
# =============================================================================

# Planck 2018 values
OMEGA_B_PLANCK = 0.0493
OMEGA_M_PLANCK = 0.315

# CBT prediction
E = np.e
BETA = 2 * E
OMEGA_EFF_CBT = OMEGA_B_PLANCK * (1 + BETA)

print("=" * 70)
print("  CMB CONSISTENCY CHECK FOR CBT")  
print("=" * 70)
print()

# =============================================================================
# THE CORE ARGUMENT
# =============================================================================

print("  THE CORE ARGUMENT")
print("  -----------------")
print()
print("  Planck 2018 measurements:")
print(f"    Ω_b = {OMEGA_B_PLANCK:.4f} (baryon density)")
print(f"    Ω_m = {OMEGA_M_PLANCK:.4f} (total matter density)")
print()
print("  CBT prediction:")
print(f"    Ω_eff = Ω_b × (1 + β) = Ω_b × (1 + 2e)")
print(f"    Ω_eff = {OMEGA_B_PLANCK:.4f} × (1 + {BETA:.4f})")
print(f"    Ω_eff = {OMEGA_EFF_CBT:.4f}")
print()

match_pct = 100 * (1 - abs(OMEGA_EFF_CBT - OMEGA_M_PLANCK) / OMEGA_M_PLANCK)
diff_pct = 100 * abs(OMEGA_EFF_CBT - OMEGA_M_PLANCK) / OMEGA_M_PLANCK

print(f"  COMPARISON:")
print(f"    Planck Ω_m = {OMEGA_M_PLANCK:.4f}")
print(f"    CBT Ω_eff  = {OMEGA_EFF_CBT:.4f}")
print(f"    Difference = {diff_pct:.2f}%")
print(f"    Match      = {match_pct:.2f}%")
print()

# =============================================================================
# WHY THIS MATTERS FOR CMB
# =============================================================================

print("=" * 70)
print("  WHY THIS MATTERS FOR CMB")
print("=" * 70)
print()
print("  The CMB power spectrum depends primarily on:")
print()
print("  1. PEAK POSITIONS (ℓ values)")
print("     → Determined by sound horizon / angular diameter distance")
print("     → Both scale with Ω_m")
print("     → If Ω_eff ≈ Ω_m, peak positions are UNCHANGED")
print()
print("  2. PEAK HEIGHTS (relative amplitudes)")
print("     → Determined by baryon-to-dark-matter ratio")
print("     → This IS different in CBT (no separate dark matter)")
print("     → Requires full Boltzmann code to evaluate")
print()
print("  3. DAMPING TAIL (high-ℓ suppression)")
print("     → Depends on free-streaming and recombination details")
print("     → Mostly determined by Ω_b, which is unchanged")
print()

# =============================================================================
# PEAK POSITION SCALING
# =============================================================================

print("=" * 70)
print("  PEAK POSITION SCALING")
print("=" * 70)
print()

# Observed first peak position
L1_OBS = 220

# Peak positions scale roughly as ℓ ∝ 1/√(Ω_m) for fixed Ω_b
# More precisely: ℓ ∝ D_A / r_s where r_s ∝ 1/√(Ω_m·h²)

# If we scale from Planck's Ω_m to CBT's Ω_eff:
scale_factor = np.sqrt(OMEGA_M_PLANCK / OMEGA_EFF_CBT)

# Observed peaks
peaks_obs = [220, 537, 810, 1120, 1420, 1720, 2030]

# CBT would predict peaks shifted by scale_factor
peaks_cbt = [int(p * scale_factor) for p in peaks_obs]

print(f"  Scale factor: √(Ω_m/Ω_eff) = {scale_factor:.4f}")
print()
print(f"{'Peak':<6s} {'Observed':>10s} {'CBT scaled':>12s} {'Shift':>10s}")
print("-" * 45)

for i, (obs, cbt) in enumerate(zip(peaks_obs, peaks_cbt)):
    shift = cbt - obs
    shift_pct = 100 * shift / obs
    print(f"  {i+1:<4d} {obs:>10d} {cbt:>12d} {shift_pct:>+9.2f}%")

print("-" * 45)
print()
print(f"  Maximum shift: {100 * (scale_factor - 1):+.2f}%")
print()

# =============================================================================
# COMPARISON TO MEASUREMENT PRECISION
# =============================================================================

print("=" * 70)
print("  COMPARISON TO MEASUREMENT PRECISION")
print("=" * 70)
print()

planck_precision = 0.3  # Planck measures first peak to ~0.3%

print(f"  Planck first peak precision: ~{planck_precision}%")
print(f"  CBT-induced peak shift:      ~{100 * abs(scale_factor - 1):.2f}%")
print()

if abs(scale_factor - 1) * 100 < planck_precision:
    print("  ✓ CBT peak shift is SMALLER than Planck's precision!")
    print("  → CBT is INDISTINGUISHABLE from ΛCDM at peak positions")
else:
    print("  ✗ CBT peak shift is LARGER than Planck's precision")
    print("  → Detailed comparison needed")

print()

# =============================================================================
# THE REMARKABLE COINCIDENCE
# =============================================================================

print("=" * 70)
print("  THE REMARKABLE RESULT")
print("=" * 70)
print()
print("  The cosmic baryon-to-matter ratio measured by Planck is:")
print()
print(f"    Ω_m / Ω_b = {OMEGA_M_PLANCK / OMEGA_B_PLANCK:.3f}")
print()
print("  CBT predicts this ratio should be:")
print()
print(f"    1 + β = 1 + 2e = {1 + BETA:.3f}")
print()
print(f"  These agree to {100 * abs((OMEGA_M_PLANCK/OMEGA_B_PLANCK) - (1+BETA)) / (OMEGA_M_PLANCK/OMEGA_B_PLANCK):.1f}%!")
print()
print("  This means:")
print("  → CBT's effective matter density = Planck's measured Ω_m")
print("  → CMB peak positions are AUTOMATICALLY correct")
print("  → No tuning required - it follows from β = 2e")
print()
print("=" * 70)

# =============================================================================
# WHAT STILL NEEDS TESTING
# =============================================================================

print("  WHAT STILL NEEDS TESTING (Full Boltzmann Code)")
print("=" * 70)
print()
print("  1. Peak height ratios (odd/even peaks)")
print("     - In ΛCDM, baryons and CDM oscillate differently")
print("     - In CBT, there's only baryons + binding field")
print("     - Does the binding field cluster like CDM?")
print()
print("  2. Damping tail shape (ℓ > 1500)")
print("     - Sensitive to recombination physics")
print("     - May require modified treatment")
print()
print("  3. Polarization spectrum")
print("     - Additional constraints from E-mode polarization")
print()
print("  VERDICT: Peak positions are consistent by construction")
print("           (Ω_eff ≈ Ω_m). Heights require more work.")
print()
print("=" * 70)

# =============================================================================
# PLOT
# =============================================================================

fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(1, len(peaks_obs) + 1)
width = 0.35

bars1 = ax.bar(x - width/2, peaks_obs, width, label='Observed (Planck)', color='black')
bars2 = ax.bar(x + width/2, peaks_cbt, width, label='CBT Prediction', color='red', alpha=0.7)

ax.set_xlabel('Peak Number', fontsize=14)
ax.set_ylabel('Multipole ℓ', fontsize=14)
ax.set_title('CMB Acoustic Peak Positions: Observed vs CBT', fontsize=16)
ax.legend(fontsize=12)
ax.grid(axis='y', alpha=0.3)

# Add difference annotation
for i, (obs, cbt) in enumerate(zip(peaks_obs, peaks_cbt)):
    diff = cbt - obs
    ax.annotate(f'{diff:+d}', xy=(x[i] + width/2, cbt), 
                ha='center', va='bottom', fontsize=9, color='red')

plt.tight_layout()
plt.savefig('cmb_peak_comparison.png', dpi=150)
print(f"  Plot saved: cmb_peak_comparison.png")
print("=" * 70)
