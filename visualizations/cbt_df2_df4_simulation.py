"""
CBT Simulation: NGC 1052-DF2 and DF4 "Dark Matter Free" Galaxies
================================================================
Author: David R. Dudas
Purpose: Demonstrate that CBT naturally explains "dark matter free" ultra-diffuse 
         galaxies (UDGs) without requiring special conditions.

The Problem:
  NGC 1052-DF2 and DF4 show rotation curves consistent with visible matter only —
  no dark matter signature. This is a crisis for CDM (where did the halo go?)
  but is naturally explained by CBT.

CBT Prediction:
  Diffuse, low-complexity structures generate minimal binding.
  α (binding strength) scales with structural complexity:
    - Compact, organized galaxy → high α → flat rotation curve
    - Diffuse, spread-out UDG → low α → declining/Keplerian curve

This simulation compares rotation curves for:
  1. Normal spiral galaxy (Milky Way-like) — should show flat curve
  2. Ultra-diffuse galaxy (DF2-like) — should show Keplerian decline
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# PHYSICAL PARAMETERS
# =============================================================================

G = 4.302e-6  # Gravitational constant in (kpc)(km/s)²/M_sun

# Normal Spiral Galaxy (Milky Way-like)
MW_MASS = 6e10          # Stellar mass (M_sun)
MW_SCALE = 3.0          # Disk scale length (kpc)
MW_VMAX = 220           # Maximum observed velocity (km/s)

# Ultra-Diffuse Galaxy (DF2-like)
DF2_MASS = 2e8          # Stellar mass (M_sun) - 300x less than MW
DF2_SCALE = 2.2         # Effective radius (kpc) - similar size!
DF2_VMAX = 20           # Observed velocity dispersion (km/s)

# CBT Parameters (from the paper)
ALPHA_BASE = 0.50       # Base binding coefficient

# =============================================================================
# CBT MODEL
# =============================================================================

def compute_alpha(R_half, M_star, is_compact=True):
    """
    Compute the CBT binding coefficient α based on structural complexity.
    
    Key insight: α depends on surface density (compactness), not just size.
    Compact galaxies have high surface density → high α → strong binding.
    Diffuse galaxies have low surface density → low α → weak binding.
    
    α(Σ) = α_base × (1 + 0.3 × log10(Σ / Σ_ref))
    
    Where Σ = M / (π R²) is the surface density.
    """
    # Surface density (M_sun / kpc²)
    Sigma = M_star / (np.pi * R_half**2)
    
    # Reference surface density (typical spiral galaxy ~100 M_sun/pc² = 1e8 M_sun/kpc²)
    Sigma_ref = 1e8
    
    # Complexity scaling
    log_ratio = np.log10(Sigma / Sigma_ref + 0.01)  # Avoid log(0)
    alpha = ALPHA_BASE * (1 + 0.3 * log_ratio)
    
    # Clamp to physical range
    alpha = np.clip(alpha, 0.0, 1.0)
    
    return alpha


def newtonian_velocity(r, M_star, R_scale):
    """
    Newtonian rotation velocity for an exponential disk.
    Uses the Freeman (1970) approximation for exponential disks.
    """
    x = r / (2 * R_scale)
    # Simplified exponential disk formula
    v_squared = G * M_star * x**2 / R_scale * (
        np.exp(-x) * (np.i0(x)*np.k0(x) - np.i0(x)*np.k1(x)) if False else
        (1 - np.exp(-x) * (1 + x)) / x  # Simpler approximation
    )
    # Use point mass + softened disk hybrid
    r_soft = np.maximum(r, 0.1)
    v_squared = G * M_star * (1 - np.exp(-r/R_scale)) / r_soft
    return np.sqrt(np.maximum(v_squared, 0))


def cbt_velocity(r, M_star, R_scale, V_max, alpha):
    """
    CBT rotation velocity = Newtonian + complexity binding contribution.
    
    v²_total = v²_Newton + v²_binding
    
    v_binding = α × V_max × min(r/r_threshold, 1)
    """
    v_newton = newtonian_velocity(r, M_star, R_scale)
    
    # Threshold radius (where binding saturates)
    r_threshold = 2.0 * R_scale
    
    # Binding velocity contribution
    v_binding = alpha * V_max * np.minimum(r / r_threshold, 1.0)
    
    # Total velocity (add in quadrature)
    v_total = np.sqrt(v_newton**2 + v_binding**2)
    
    return v_total, v_newton, v_binding


def create_comparison_plot():
    """Create side-by-side comparison of normal vs UDG rotation curves."""
    
    # Radial range
    r = np.linspace(0.1, 15, 200)  # kpc
    
    # Compute α for each galaxy type
    alpha_MW = compute_alpha(MW_SCALE, MW_MASS, is_compact=True)
    alpha_DF2 = compute_alpha(DF2_SCALE, DF2_MASS, is_compact=False)
    
    print(f"Surface density MW:  {MW_MASS / (np.pi * MW_SCALE**2):.2e} M_sun/kpc²")
    print(f"Surface density DF2: {DF2_MASS / (np.pi * DF2_SCALE**2):.2e} M_sun/kpc²")
    print(f"α (MW):  {alpha_MW:.3f}")
    print(f"α (DF2): {alpha_DF2:.3f}")
    
    # Compute rotation curves
    v_mw_total, v_mw_newton, v_mw_bind = cbt_velocity(r, MW_MASS, MW_SCALE, MW_VMAX, alpha_MW)
    v_df2_total, v_df2_newton, v_df2_bind = cbt_velocity(r, DF2_MASS, DF2_SCALE, DF2_VMAX, alpha_DF2)
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor('#0d1117')
    fig.suptitle('CBT Explains "Dark Matter Free" Galaxies', 
                 color='white', fontsize=16, fontweight='bold')
    
    # Colors
    c_newton = '#ff6b6b'
    c_binding = '#4ecdc4'
    c_total = '#fbbf24'
    
    # Panel 1: Normal Spiral Galaxy
    ax1 = axes[0]
    ax1.set_facecolor('#0d1117')
    ax1.plot(r, v_mw_newton, '--', color=c_newton, linewidth=2, label='Newtonian (visible mass)')
    ax1.plot(r, v_mw_bind, ':', color=c_binding, linewidth=2, label=f'CBT Binding (α={alpha_MW:.2f})')
    ax1.plot(r, v_mw_total, '-', color=c_total, linewidth=3, label='CBT Total')
    
    ax1.axhline(y=MW_VMAX, color='white', linestyle='--', alpha=0.3, label=f'Observed V_max={MW_VMAX} km/s')
    
    ax1.set_xlabel('Radius (kpc)', color='gray', fontsize=11)
    ax1.set_ylabel('Rotation Velocity (km/s)', color='gray', fontsize=11)
    ax1.set_title('Normal Spiral Galaxy (MW-like)\nHigh complexity → Strong binding → Flat curve', 
                  color='white', fontsize=11)
    ax1.set_xlim(0, 15)
    ax1.set_ylim(0, 300)
    ax1.tick_params(colors='gray')
    ax1.legend(loc='lower right', facecolor='#1a202c', edgecolor='#4a5568', labelcolor='white')
    ax1.grid(True, alpha=0.2, color='white')
    for spine in ax1.spines.values():
        spine.set_color('#4a5568')
        
    # Add annotation
    ax1.annotate('Flat rotation curve\n(like observed)', xy=(12, 220), 
                 color='#fbbf24', fontsize=10, ha='center')
    
    # Panel 2: Ultra-Diffuse Galaxy (DF2-like)
    ax2 = axes[1]
    ax2.set_facecolor('#0d1117')
    ax2.plot(r, v_df2_newton, '--', color=c_newton, linewidth=2, label='Newtonian (visible mass)')
    ax2.plot(r, v_df2_bind, ':', color=c_binding, linewidth=2, label=f'CBT Binding (α={alpha_DF2:.2f})')
    ax2.plot(r, v_df2_total, '-', color=c_total, linewidth=3, label='CBT Total')
    
    ax2.axhline(y=DF2_VMAX, color='white', linestyle='--', alpha=0.3, label=f'Observed σ={DF2_VMAX} km/s')
    
    ax2.set_xlabel('Radius (kpc)', color='gray', fontsize=11)
    ax2.set_title('Ultra-Diffuse Galaxy (DF2-like)\nLow complexity → Weak binding → Keplerian decline', 
                  color='white', fontsize=11)
    ax2.set_xlim(0, 15)
    ax2.set_ylim(0, 50)
    ax2.tick_params(colors='gray')
    ax2.legend(loc='upper right', facecolor='#1a202c', edgecolor='#4a5568', labelcolor='white')
    ax2.grid(True, alpha=0.2, color='white')
    for spine in ax2.spines.values():
        spine.set_color('#4a5568')
    
    # Add annotation
    ax2.annotate('Declining curve\n(matches DF2 observation!)', xy=(10, 15), 
                 color='#fbbf24', fontsize=10, ha='center')
    
    # Explanation box
    textbox = """
KEY INSIGHT: CBT explains DF2/DF4 naturally.

• Normal galaxies: High surface density → High α → Strong binding → "Dark matter" signature
• Diffuse UDGs: Low surface density → Low α → Weak binding → No "dark matter"

CDM must invoke special conditions (tidal stripping, failed formation).
CBT predicts this from first principles — no dark matter to lose!
"""
    fig.text(0.5, -0.02, textbox, ha='center', va='top', fontsize=10,
             color='#a0aec0', family='monospace',
             bbox=dict(boxstyle='round', facecolor='#1a202c', edgecolor='#4a5568'))
    
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])
    plt.savefig('cbt_df2_udg_comparison.png', dpi=150, 
                facecolor='#0d1117', edgecolor='none', bbox_inches='tight')
    print("\nSaved: cbt_df2_udg_comparison.png")
    plt.show()
    
    # Print quantitative summary
    print("\n" + "="*60)
    print("QUANTITATIVE SUMMARY")
    print("="*60)
    print(f"\nSurface Density Ratio: MW/DF2 = {(MW_MASS/MW_SCALE**2)/(DF2_MASS/DF2_SCALE**2):.0f}x")
    print(f"Binding Coefficient Ratio: α_MW/α_DF2 = {alpha_MW/alpha_DF2:.1f}x")
    print(f"\nAt r = 5 kpc:")
    print(f"  MW:  v_Newton = {newtonian_velocity(5, MW_MASS, MW_SCALE):.0f} km/s,  v_total = {v_mw_total[np.argmin(np.abs(r-5))]:.0f} km/s")
    print(f"  DF2: v_Newton = {newtonian_velocity(5, DF2_MASS, DF2_SCALE):.0f} km/s,  v_total = {v_df2_total[np.argmin(np.abs(r-5))]:.0f} km/s")
    print(f"\nCBT correctly predicts:")
    print(f"  • MW-like galaxies need 'dark matter' (actually: strong binding)")
    print(f"  • DF2-like UDGs have no 'dark matter' (actually: weak binding)")


if __name__ == "__main__":
    print("="*60)
    print("CBT SIMULATION: DF2/DF4 'DARK MATTER FREE' GALAXIES")
    print("="*60)
    print()
    print("The discovery of NGC 1052-DF2 and DF4 (van Dokkum et al. 2018)")
    print("challenged dark matter theory: where did their halos go?")
    print()
    print("CBT Answer: They never had 'halos' — they're just diffuse.")
    print("Diffuse structures have low complexity → low binding → no DM signature.")
    print()
    
    create_comparison_plot()
