import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import re

def load_all_csvs():
    """Load multi-follicle CSV files"""
    files = sorted(glob.glob("output/phi_t*.csv"))
    if not files:
        raise FileNotFoundError("No CSV files found in 'output/'")
    
    times = []
    m_grid = None
    phi_f0_all, phi_f1_all, phi_f2_all = [], [], []
    
    for f in files:
        df = pd.read_csv(f)
        
        # Extract time from filename
        step_match = re.search(r'phi_t(\d+)\.csv$', f)
        step = int(step_match.group(1)) if step_match else 0
        time = step * 0.012  # CFL dt ≈ 0.016
        
        times.append(time)
        phi_f0_all.append(df['phi_f0'].values)
        phi_f1_all.append(df['phi_f1'].values) 
        phi_f2_all.append(df['phi_f2'].values)
        
        if m_grid is None:
            m_grid = df['m'].values
            
        print(f"Loaded {f}: t={time:.3f}, Nf=3")
    
    return np.array(times), m_grid, np.array(phi_f0_all), np.array(phi_f1_all), np.array(phi_f2_all)

def plot_evolution():
    """Plot all 3 follicles over time"""
    times, m, phi_f0, phi_f1, phi_f2 = load_all_csvs()
    
    plt.figure(figsize=(12, 8))
    n_plots = min(10, len(times))
    indices = np.linspace(0, len(times)-1, n_plots, dtype=int)
    
    for i in indices:
        plt.plot(m, phi_f0[i], 'b-', alpha=0.8, linewidth=2, label=f't={times[i]:.3f}' if i==0 else "")
        plt.plot(m, phi_f1[i], 'r-', alpha=0.8, linewidth=2)
        plt.plot(m, phi_f2[i], 'g-', alpha=0.8, linewidth=2)
    
    plt.xlabel('Maturity $m$', fontsize=14)
    plt.ylabel('Density $\\phi(m)$', fontsize=14)
    plt.title('Multi-Follicle Evolution (Nf=3)', fontsize=16)
    plt.legend(['Follicle 0', 'Follicle 1', 'Follicle 2'])
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('multi_follicle_evolution.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_individual_follicles():
    """Separate plot for each follicle"""
    times, m, phi_f0, phi_f1, phi_f2 = load_all_csvs()
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # Follicle 0
    n_times = min(8, len(times))
    indices = np.linspace(0, len(times)-1, n_times, dtype=int)
    for i in indices:
        axes[0].plot(m, phi_f0[i], 'b-', alpha=0.7)
    axes[0].set_title('Follicle 0')
    axes[0].set_xlabel('Maturity $m$')
    axes[0].set_ylabel('$\\phi(m)$')
    axes[0].grid(True, alpha=0.3)
    
    # Follicle 1  
    for i in indices:
        axes[1].plot(m, phi_f1[i], 'r-', alpha=0.7)
    axes[1].set_title('Follicle 1')
    axes[1].set_xlabel('Maturity $m$')
    axes[1].grid(True, alpha=0.3)
    
    # Follicle 2
    for i in indices:
        axes[2].plot(m, phi_f2[i], 'g-', alpha=0.7)
    axes[2].set_title('Follicle 2')
    axes[2].set_xlabel('Maturity $m$')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('individual_follicles.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_mass_evolution():
    """Mass evolution for each follicle"""
    times, m, phi_f0, phi_f1, phi_f2 = load_all_csvs()
    
    mass_f0 = [np.trapezoid(phi, m) for phi in phi_f0]
    mass_f1 = [np.trapezoid(phi, m) for phi in phi_f1]
    mass_f2 = [np.trapezoid(phi, m) for phi in phi_f2]
    
    plt.figure(figsize=(10, 6))
    plt.plot(times, mass_f0, 'bo', linewidth=3, label='Follicle 0', markersize=6)
    plt.plot(times, mass_f1, 'rs', linewidth=3, label='Follicle 1', markersize=6)
    plt.plot(times, mass_f2, 'g^', linewidth=3, label='Follicle 2', markersize=6)
    
    plt.xlabel('Time $t$', fontsize=14)
    plt.ylabel('Mass per follicle', fontsize=14)
    plt.title('Individual Follicle Mass Evolution', fontsize=16)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('follicle_masses.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    print("Loading multi-follicle CSV files (m,phi_f0,phi_f1,phi_f2,g,lambda,p)...")
    plot_evolution()
    plot_individual_follicles()
    plot_mass_evolution()
    print("Plots saved:")
    print("- multi_follicle_evolution.png")
    print("- individual_follicles.png") 
    print("- follicle_masses.png")

