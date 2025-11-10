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
frames = phi_files[:100]

fig, ax = plt.subplots()
# Initialize with first frame
phi_1d = np.loadtxt(frames[0])
phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))
im = ax.imshow(phi.T, extent=[coords_x[0], coords_x[-1], coords_y[0], coords_y[-1]],
                         origin='lower', cmap='bwr', animated=True)
cbar = plt.colorbar(im, ax=ax, label='phi')

def update(frame_idx):
    filename = frames[frame_idx]
    phi_1d = np.loadtxt(filename)
    phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))
    # Compute new vmin, vmax each frame
    vmin = np.min(phi)
    vmax = np.max(phi)
    # Update image data
    im.set_data(phi.T)
    # Update normalization to reflect current vmin/vmax
    im.set_norm(colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=vmin, vmax=vmax, base=10))
    # Optionally update colorbar if needed (not shown here, but possible)
    return [im]

ani = FuncAnimation(fig, update, frames=len(frames), interval=100, blit=False)
plt.show()

