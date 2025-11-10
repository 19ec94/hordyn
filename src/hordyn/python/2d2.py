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
step_count_end = 100
frames = phi_files[step_count_start:step_count_end]

fig, ax = plt.subplots(figsize=(8, 1))
pcm = None  # QuadMesh placeholder
cbar = None # Colorbar placeholder

def update(frame_idx):
    global pcm, cbar
    filepath = frames[frame_idx]
    phi_1d = np.loadtxt(filepath)
    phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))

    vmin = np.min(phi[np.nonzero(phi)]) if np.any(phi > 0) else 1e-10
    vmax = np.max(phi)

    if pcm is None:
        norm = colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=vmin, vmax=vmax, base=10)
        pcm = ax.pcolormesh(X, Y, phi, shading='auto', cmap='bwr', norm=norm)
        cbar = fig.colorbar(pcm, ax=ax, label='phi')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    else:
        norm = colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=vmin, vmax=vmax, base=10)
        pcm.set_norm(norm)
        pcm.set_array(phi.ravel())
        cbar.update_normal(pcm)

    ax.set_title(f'Phi heatmap from {os.path.basename(filepath)}')
    return pcm, cbar

ani = FuncAnimation(fig, update, frames=len(frames), interval=100, blit=False, repeat=False)
plt.show()

