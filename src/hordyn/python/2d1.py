import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import glob
import os
import re
from matplotlib.animation import FuncAnimation

folder = 'output/f0'
file_pattern = os.path.join(folder, 'phi*.txt')

def numerical_sort(value):
    match = re.search(r'phi(\d+)\.txt', value)
    if match:
        return int(match.group(1))
    else:
        return -1

phi_files = sorted(glob.glob(file_pattern), key=numerical_sort)

coords_x = np.loadtxt(os.path.join(folder, 'coords_x.txt'))
coords_y = np.loadtxt(os.path.join(folder, 'coords_y.txt'))
X, Y = np.meshgrid(coords_x, coords_y)

step_count_start = 0
step_count_end = 600
frames = phi_files[step_count_start:step_count_end]

fig, ax = plt.subplots(figsize=(8,1))
pcm = None  # QuadMesh object placeholder

# Pre-read min and max of phi across all selected frames to fix color scale
all_min = float('inf')
all_max = float('-inf')
for filepath in frames:
    phi_1d = np.loadtxt(filepath)
    phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))
    all_min = min(all_min, np.min(phi))
    all_max = max(all_max, np.max(phi))

# Create initial plot with first frame
phi_1d_init = np.loadtxt(frames[0])
phi_init = phi_1d_init.reshape((len(coords_y)-1, len(coords_x)-1))
pcm = ax.pcolormesh(X, Y, phi_init, shading='auto', cmap='bwr',
                    norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=all_min, vmax=all_max, base=10))
cbar = fig.colorbar(pcm, ax=ax, label='phi')
ax.set_xlabel('x')
ax.set_ylabel('y')
title_text = ax.set_title(f'Phi heatmap from {os.path.basename(frames[0])}')

def update(frame_idx):
    filepath = frames[frame_idx]
    phi_1d = np.loadtxt(filepath)
    phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))
    # Update QuadMesh data
    pcm.set_array(phi.ravel())
    # Color limits fixed globally - no need to reset set_clim here
    title_text.set_text(f'Phi heatmap from {os.path.basename(filepath)}')
    return pcm, title_text

ani = FuncAnimation(fig, update, frames=len(frames), interval=100, blit=False, repeat=False)

plt.show()

