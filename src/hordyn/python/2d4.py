import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import glob
import os
import re
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm


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

step_count_start = 0 
step_count_end = 1500
frames = phi_files[step_count_start:step_count_end]

fig, ax = plt.subplots(figsize=(1, 1))

# Load first frame
phi_1d = np.loadtxt(frames[0])
phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))

vmin, vmax = np.min(phi), np.max(phi)
boundaries = [vmin, vmin + 0.2*(vmax-vmin), vmin + 0.70*(vmax-vmin), vmax]
# Create your custom segmented colormap
#colors = [
#    (0.0, 'blue'),
#    #(0.3, 'white'),
#    (0.5, 'green'),
#    #(0.6, 'yellow'),
#    (1.0, 'red')
#]
colors = ['blue', 'green', 'red']
palette = ['#000fff', '#90ff70', '#ee0000']
segmented_cmap = LinearSegmentedColormap.from_list('custom_palette', palette, N=256)
norm = BoundaryNorm(boundaries, ncolors=segmented_cmap.N, clip=True)


# Create imshow plot with initial linear normalization
im = ax.imshow(phi.T,
               extent=[coords_x[0], coords_x[-1], coords_y[0], coords_y[-1]],
               origin='lower',
               cmap=segmented_cmap,
               #vmin=phi.min(),
               #vmax=phi.max(),
               norm = norm,
               aspect='auto')  # Aspect set for rectangular visualization

cbar = fig.colorbar(im, ax=ax, label='phi')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal')
title_text = ax.set_title(f'Phi heatmap from {os.path.basename(frames[0])}')

# Code added - start
Ys = 0.3  # Horizontal line y position


# Draw horizontal line spanning whole x axis at y=Ys
hline = ax.axhline(y=Ys, color='black', linewidth=2)

# Draw vertical lines at specified x positions, only spanning from y=0 to y=Ys
# ymin and ymax are relative positions in axis coordinates (0 to 1)
#v_positions = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
v_positions = [0.5, 1.0]
colors = ['gray', 'black'] * (len(v_positions)//2 + 1)
widths = [1, 2] * (len(v_positions)//2 + 1)

v_lines = []
for x, c, w in zip(v_positions, colors, widths):
    line = ax.axvline(x=x, ymin=0.0, ymax=Ys, color=c, linewidth=w)
    v_lines.append(line)
# Code added - end
def update(frame_idx):
    filename = frames[frame_idx]
    phi_1d = np.loadtxt(filename)
    phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))

    vmin = np.min(phi)
    vmax = np.max(phi)

    im.set_data(phi)
    im.set_clim(vmin, vmax)  # Dynamic linear color limits
    title_text.set_text(f'Phi heatmap from {os.path.basename(filename)}')

    return [im, title_text]

ani = FuncAnimation(fig, update, frames=len(frames), interval=100, blit=False, repeat=False)

plt.show()

