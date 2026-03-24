import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import PowerNorm

def plot_phi_heatmap(cmap_name='inferno'):
    files = ["pcos/phi0.csv", "pcos/phi1.csv", "pcos/phi2.csv", "pcos/phi3.csv"]
    epsilon = 1e-10
    
    # Load real time values from case_log.csv
    case_log = pd.read_csv("pcos/case_log.csv")
    real_times = case_log['t'].values  # Shape: (n_timesteps,)
    
    # Global bounds computation
    all_matrices = []
    for fname in files:
        df = pd.read_csv(fname)
        phi_cols = [col for col in df.columns if col.startswith("phi_t")]
        phi_matrix = df[phi_cols].values + epsilon
        all_matrices.append(phi_matrix)
    
    global_min = epsilon
    global_max = max(np.max(mat) for mat in all_matrices)
    shared_norm = PowerNorm(gamma=0.3, vmin=global_min, vmax=global_max)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.ravel()
    
    for ax, fname in zip(axes, files):
        df = pd.read_csv(fname)
        phi_cols = [col for col in df.columns if col.startswith("phi_t")]
        m_values = df['m'].values
        phi_matrix = df[phi_cols].values + epsilon
        
        # Use real time values: [t_min, t_max, m_min, m_max]
        extent = [real_times[0], real_times[-1], m_values[0], m_values[-1]]
        
        im = ax.imshow(phi_matrix, aspect='auto', origin='lower', 
                       cmap=cmap_name,
                       extent=extent,  # Real t values on x-axis!
                       norm=shared_norm)
        
        ax.set_xlabel('t (real time, days)')
        ax.set_ylabel('m')
        phi_index = int(fname[8]) if len(fname) > 3 else 0
        ax.set_title(f'$\\phi_{{{phi_index}}}$')
        plt.colorbar(im, ax=ax, shrink=0.8, label='φ (power scale)')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_phi_heatmap('gist_heat_r')
    # plot_phi_heatmap('RdYlBu_r')

