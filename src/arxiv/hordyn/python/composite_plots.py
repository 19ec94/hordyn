import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import glob, os, re
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec

# ---------- Load Data ----------

folder_anim = 'output/fullmodel/2012_tau1.0_v2/f0'
file_pattern = os.path.join(folder_anim, 'phi*.txt')
summary_file_anim = os.path.join(folder_anim, 'summary.csv')
coords_x = np.loadtxt(os.path.join(folder_anim, 'coords_x.txt'))
coords_y = np.loadtxt(os.path.join(folder_anim, 'coords_y.txt'))

# Helper function for sorting frames numerically
def numerical_sort(value):
    match = re.search(r'phi(\d+)\.txt', value)
    return int(match.group(1)) if match else -1

phi_files = sorted(glob.glob(file_pattern), key=numerical_sort)
df_anim = pd.read_csv(summary_file_anim)
time_stamps = df_anim['time']

# For time series plots
data_path_1 = 'output/fullmodel/2012_tau1.0_v2/f0/summary.csv'
data_path_2 = 'output/fullmodel/2012_tau1.0/f0/summary.csv'
df_1 = pd.read_csv(data_path_1)
df_2 = pd.read_csv(data_path_2)

# ---------- Define Plot Variables ----------

plot_vars = [
    {'column': 'phi_max', 'ylabel': 'phi_max', 'title': 'Time vs phi_max', 'color': None},
    {'column': 'mass', 'ylabel': 'Mass', 'title': 'Time vs Mass', 'color': 'green'},
    # add more if needed...
]

# --------- Set up Combined Figure ---------

fig = plt.figure(figsize=(16, 12))
gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1])

# Top (Animation) - spans both columns
ax_anim = fig.add_subplot(gs[0, :])

# Bottom (Static time-series, side-by-side)
ax_ts1 = fig.add_subplot(gs[1, 0])
ax_ts2 = fig.add_subplot(gs[1, 1])

# --------- Prepare Initial Frame for Animation ---------

phi_1d = np.loadtxt(phi_files[0])
phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))
palette = ['#000fff', '#90ff70', '#ee0000']
color_steps = 25

def get_boundaries(data, steps=10):
    vmin, vmax = data.min(), data.max()
    blue_end = vmin + 0.1 * (vmax - vmin)
    red_start = vmin + 0.95 * (vmax - vmin)
    boundaries = np.linspace(blue_end, red_start, steps)
    boundaries = np.concatenate(([vmin, blue_end], boundaries, [red_start, vmax]))
    return boundaries

segmented_cmap = LinearSegmentedColormap.from_list('red_green_blue', palette, N=color_steps)
boundaries = get_boundaries(phi, color_steps)
num_color_intervals = len(boundaries) -1
norm = BoundaryNorm(boundaries, ncolors=num_color_intervals, clip=True)

im = ax_anim.imshow(phi,
    extent=[coords_x[0], coords_x[-1], coords_y[0], coords_y[-1]],
    origin='lower',
    cmap=segmented_cmap,
    norm=norm,
    aspect='auto')

cbar = fig.colorbar(im, ax=ax_anim, label='phi')
ax_anim.set_xlabel('x')
ax_anim.set_ylabel('y')
ax_anim.set_aspect('equal')
title_text = ax_anim.text(0.5, 1.02, f'Phi heatmap at time = {time_stamps.iloc[0]}', transform=ax_anim.transAxes, ha='center')
# Horizontal and vertical lines:
Ys = 0.3
ax_anim.axhline(y=Ys, color='black', linewidth=2)
v_positions = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]
colors_lines = ['gray', 'black'] * (len(v_positions)//2 + 1)
widths = [1, 2] * (len(v_positions)//2 + 1)
for x, c, w in zip(v_positions, colors_lines, widths):
    ax_anim.axvline(x=x, ymin=0.0, ymax=Ys, color=c, linewidth=w)

# ---------- Animation Update ----------

def update(frame_idx):
    filename = phi_files[frame_idx]
    phi_1d = np.loadtxt(filename)
    phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))
    boundaries = get_boundaries(phi, color_steps)
    num_color_intervals = len(boundaries) - 1
    norm = BoundaryNorm(boundaries, ncolors=num_color_intervals, clip=True)
    im.set_norm(norm)
    im.set_data(phi)
    title_text.set_text(f'Phi heatmap at time = {time_stamps.iloc[frame_idx]}')
    return [im, title_text]

ani = FuncAnimation(fig, update, frames=len(phi_files), interval=100, blit=False, repeat=False)

# ---------- Time-Series Plots (Below) ----------

v1 = plot_vars[0]
ax = ax_ts1
ax.plot(df_1['time'], df_1[v1['column']], label="tau=1.0")
ax.plot(df_2['time'], df_2[v1['column']], label="tau=1.2")
ax.set_xlabel('Time')
ax.set_ylabel(v1['ylabel'])
ax.set_title(v1['title'])
ax.grid(True)
ax.legend()

v2 = plot_vars[1]
ax = ax_ts2
ax.plot(df_1['time'], df_1[v2['column']], label="tau=1.0")
ax.plot(df_2['time'], df_2[v2['column']], label="tau=1.2")
ax.set_xlabel('Time')
ax.set_ylabel(v2['ylabel'])
ax.set_title(v2['title'])
ax.grid(True)
ax.legend()

plt.tight_layout()
plt.show()

