import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def main():
    files = ["pcos/phi0.csv", "pcos/phi1.csv", "pcos/phi2.csv", "pcos/phi3.csv"]
    x_col = "m"
    colors = ['blue', 'orange', 'green', 'red']
    
    # Load real time values from case_log.csv
    case_log = pd.read_csv("pcos/case_log.csv")
    real_times = case_log['t'].values
    
    # Read all files and get phi columns
    dfs = [pd.read_csv(f) for f in files]
    phi_col_lists = [[col for col in df.columns if col.startswith("phi_t")] for df in dfs]
    
    # Use SHORTEST length across ALL data sources (no assert needed)
    total_frames = min(len(real_times), *[len(phi_cols) for phi_cols in phi_col_lists])
    print(f"Using {total_frames} frames (case_log: {len(real_times)}, phi files: {len(phi_cols) for phi_cols in phi_col_lists})")
    
    # Create single plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Initialize one line per file
    lines = []
    for i, (df, phi_cols, color) in enumerate(zip(dfs, phi_col_lists, colors)):
        line, = ax.plot(df[x_col], df[phi_cols[0]], label=f"$\\phi_{{{i}}}$", 
                       color=color, linewidth=2)
        lines.append(line)
    
    ax.set_xlabel(x_col)
    ax.set_ylabel("φ")
    ax.set_title("φ evolution - all files")
    ax.grid(True)
    ax.legend()
    
    def animate(frame):
        # Update all lines to current time slice
        for i, (line, df, phi_cols) in enumerate(zip(lines, dfs, phi_col_lists)):
            line.set_data(df[x_col], df[phi_cols[frame]])
        
        # Dynamic title with REAL TIME from case_log.csv
        t_val = real_times[frame]
        ax.set_title(f"φ evolution at t={t_val:.6f}")
        
        # Auto-scale
        ax.relim()
        ax.autoscale_view()
        
        return lines + [ax.title]
    
    # Create animation with safe frame count
    anim = FuncAnimation(fig, animate, frames=total_frames, interval=100, 
                         blit=False, repeat=True)
    
    plt.tight_layout()
    plt.show()
    return anim

if __name__ == "__main__":
    main()

