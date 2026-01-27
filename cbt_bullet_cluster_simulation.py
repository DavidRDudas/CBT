"""
Complexity Binding Theory (CBT) Simulation: Bullet Cluster Collision
=====================================================================
Author: David R. Dudas
Purpose: Demonstrate how a "binding force" (emergent from structural complexity)
         produces galaxy separation during cluster collisions, similar to what
         is observed in the Bullet Cluster.

Hypothesis: Standard Dark Matter assumes gravity from invisible mass.
            CBT proposes that effective gravity emerges from structural complexity
            (low entropy). Organized systems generate "binding energy" to maintain
            coherence, while disorganized systems (gas) do not.

This simulation runs two scenarios side-by-side:
  A) Newtonian: Standard 1/r² gravity only → chaotic merger
  B) CBT: Standard gravity + binding force → clusters pass through and reform
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

N_PARTICLES = 120          # Particles per cluster
CLUSTER_RADIUS = 1.2       # Initial spread
G = 0.4                    # Gravitational constant
SOFTENING = 0.4            # Softening to prevent singularities
DT = 0.025                 # Time step
N_STEPS = 500              # Total frames

# Initial conditions - head-on collision
CLUSTER_1_POS = np.array([-4.0, 0.3])
CLUSTER_2_POS = np.array([4.0, -0.3])
CLUSTER_1_VEL = np.array([1.2, 0.0])
CLUSTER_2_VEL = np.array([-1.2, 0.0])

# CBT Binding Force constant
K_BIND = 0.6

# =============================================================================
# SIMULATION CLASS
# =============================================================================

class ClusterSimulation:
    """N-body simulation with optional CBT binding force."""
    
    def __init__(self, use_binding=False):
        self.use_binding = use_binding
        self.n = N_PARTICLES
        
        # Initialize cluster 1
        theta1 = np.random.uniform(0, 2*np.pi, self.n)
        r1 = CLUSTER_RADIUS * np.sqrt(np.random.uniform(0, 1, self.n))
        pos1 = np.column_stack([
            CLUSTER_1_POS[0] + r1 * np.cos(theta1),
            CLUSTER_1_POS[1] + r1 * np.sin(theta1)
        ])
        vel1 = np.column_stack([
            CLUSTER_1_VEL[0] + np.random.normal(0, 0.08, self.n),
            CLUSTER_1_VEL[1] + np.random.normal(0, 0.08, self.n)
        ])
        
        # Initialize cluster 2
        theta2 = np.random.uniform(0, 2*np.pi, self.n)
        r2 = CLUSTER_RADIUS * np.sqrt(np.random.uniform(0, 1, self.n))
        pos2 = np.column_stack([
            CLUSTER_2_POS[0] + r2 * np.cos(theta2),
            CLUSTER_2_POS[1] + r2 * np.sin(theta2)
        ])
        vel2 = np.column_stack([
            CLUSTER_2_VEL[0] + np.random.normal(0, 0.08, self.n),
            CLUSTER_2_VEL[1] + np.random.normal(0, 0.08, self.n)
        ])
        
        # Combine
        self.positions = np.vstack([pos1, pos2])
        self.velocities = np.vstack([vel1, vel2])
        self.cluster_ids = np.array([0]*self.n + [1]*self.n)
        self.total = 2 * self.n
    
    def compute_acceleration(self):
        """Compute gravitational + binding acceleration."""
        acc = np.zeros_like(self.positions)
        
        # Gravitational acceleration (vectorized)
        for i in range(self.total):
            dx = self.positions - self.positions[i]
            r2 = np.sum(dx**2, axis=1) + SOFTENING**2
            mask = np.arange(self.total) != i
            a_grav = G * dx[mask] / (r2[mask, np.newaxis] ** 1.5)
            acc[i] = np.sum(a_grav, axis=0)
        
        # Binding force (CBT)
        if self.use_binding:
            com0 = np.mean(self.positions[self.cluster_ids == 0], axis=0)
            com1 = np.mean(self.positions[self.cluster_ids == 1], axis=0)
            
            for i in range(self.total):
                com = com0 if self.cluster_ids[i] == 0 else com1
                acc[i] -= K_BIND * (self.positions[i] - com)
        
        return acc
    
    def step(self):
        """Velocity Verlet integration step."""
        acc = self.compute_acceleration()
        self.positions += self.velocities * DT + 0.5 * acc * DT**2
        acc_new = self.compute_acceleration()
        self.velocities += 0.5 * (acc + acc_new) * DT
    
    def get_com(self, cluster_id):
        """Get center of mass for a cluster."""
        return np.mean(self.positions[self.cluster_ids == cluster_id], axis=0)
    
    def get_spread(self, cluster_id):
        """Get average distance from COM (compactness metric)."""
        mask = self.cluster_ids == cluster_id
        com = np.mean(self.positions[mask], axis=0)
        return np.mean(np.linalg.norm(self.positions[mask] - com, axis=1))

# =============================================================================
# ANIMATION
# =============================================================================

def run_simulation():
    """Run side-by-side animation comparing Newtonian vs CBT."""
    
    # Use same random seed for fair comparison
    np.random.seed(42)
    sim_newton = ClusterSimulation(use_binding=False)
    np.random.seed(42)
    sim_cbt = ClusterSimulation(use_binding=True)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor('#0d1117')
    fig.suptitle('CBT Bullet Cluster Simulation: Order vs Chaos', 
                 color='white', fontsize=16, fontweight='bold')
    
    # Colors
    c1, c2 = '#ff6b6b', '#4ecdc4'
    
    for ax, title in [(ax1, 'Scenario A: Newtonian (No Binding)'),
                      (ax2, 'Scenario B: CBT (With Binding)')]:
        ax.set_facecolor('#0d1117')
        ax.set_xlim(-12, 12)
        ax.set_ylim(-8, 8)
        ax.set_aspect('equal')
        ax.set_title(title, color='white', fontsize=12, pad=10)
        ax.tick_params(colors='#3d4450')
        for spine in ax.spines.values():
            spine.set_color('#3d4450')
        ax.set_xlabel('x', color='#666')
        ax.set_ylabel('y', color='#666')
    
    # Scatter plots
    scat1_n = ax1.scatter([], [], c=c1, s=12, alpha=0.8, label='Cluster 1')
    scat2_n = ax1.scatter([], [], c=c2, s=12, alpha=0.8, label='Cluster 2')
    scat1_c = ax2.scatter([], [], c=c1, s=12, alpha=0.8, label='Cluster 1')
    scat2_c = ax2.scatter([], [], c=c2, s=12, alpha=0.8, label='Cluster 2')
    
    # COM markers
    com1_n, = ax1.plot([], [], 'o', color='white', ms=8, mec=c1, mew=2)
    com2_n, = ax1.plot([], [], 'o', color='white', ms=8, mec=c2, mew=2)
    com1_c, = ax2.plot([], [], 'o', color='white', ms=8, mec=c1, mew=2)
    com2_c, = ax2.plot([], [], 'o', color='white', ms=8, mec=c2, mew=2)
    
    # Status text
    info_n = ax1.text(0.02, 0.98, '', transform=ax1.transAxes, color='white',
                      fontsize=9, va='top', family='monospace')
    info_c = ax2.text(0.02, 0.98, '', transform=ax2.transAxes, color='white',
                      fontsize=9, va='top', family='monospace')
    
    # Annotations
    ax1.text(0.98, 0.02, 'Clusters merge chaotically', transform=ax1.transAxes,
             color='#ff6b6b', fontsize=9, ha='right', va='bottom', style='italic')
    ax2.text(0.98, 0.02, 'Clusters pass through & reform', transform=ax2.transAxes,
             color='#4ecdc4', fontsize=9, ha='right', va='bottom', style='italic')
    
    def update(frame):
        # Step physics
        sim_newton.step()
        sim_cbt.step()
        
        # Update scatter positions
        scat1_n.set_offsets(sim_newton.positions[:N_PARTICLES])
        scat2_n.set_offsets(sim_newton.positions[N_PARTICLES:])
        scat1_c.set_offsets(sim_cbt.positions[:N_PARTICLES])
        scat2_c.set_offsets(sim_cbt.positions[N_PARTICLES:])
        
        # Update COM markers
        cn1 = sim_newton.get_com(0)
        cn2 = sim_newton.get_com(1)
        cc1 = sim_cbt.get_com(0)
        cc2 = sim_cbt.get_com(1)
        
        com1_n.set_data([cn1[0]], [cn1[1]])
        com2_n.set_data([cn2[0]], [cn2[1]])
        com1_c.set_data([cc1[0]], [cc1[1]])
        com2_c.set_data([cc2[0]], [cc2[1]])
        
        # Update info
        sep_n = np.linalg.norm(cn2 - cn1)
        sep_c = np.linalg.norm(cc2 - cc1)
        spread_n = (sim_newton.get_spread(0) + sim_newton.get_spread(1)) / 2
        spread_c = (sim_cbt.get_spread(0) + sim_cbt.get_spread(1)) / 2
        
        info_n.set_text(f'Frame: {frame}\nSeparation: {sep_n:.1f}\nSpread: {spread_n:.2f}')
        info_c.set_text(f'Frame: {frame}\nSeparation: {sep_c:.1f}\nSpread: {spread_c:.2f}')
        
        return scat1_n, scat2_n, scat1_c, scat2_c, com1_n, com2_n, com1_c, com2_c
    
    anim = FuncAnimation(fig, update, frames=N_STEPS, interval=30, blit=False)
    
    plt.tight_layout()
    plt.show()
    
    return anim

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("COMPLEXITY BINDING THEORY (CBT) SIMULATION")
    print("Bullet Cluster Collision Demonstration")
    print("=" * 60)
    print()
    print("Scenario A (Left):  Standard 1/r² gravity only")
    print("Scenario B (Right): Gravity + CBT binding force")
    print()
    print("The binding force: a = -k(r - r_COM)")
    print("This represents structural coherence from complexity.")
    print()
    print("Watch how CBT clusters pass through and maintain")
    print("their identity, while Newtonian clusters merge.")
    print("-" * 60)
    
    run_simulation()
