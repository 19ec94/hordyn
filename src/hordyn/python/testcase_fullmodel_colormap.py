import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import glob
import os
import re
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm

# Folder and file pattern for loading data files
folder = 'output/fullmodel/2012_tau1.0_v2/f0'
file_pattern = os.path.join(folder, 'phi*.txt')

# Read the timestamps
import pandas as pd
df = pd.read_csv(os.path.join(folder, 'summary.csv'))
time_stamps = df['time'] 

def numerical_sort(value):
    """
    Sort filenames numerically based on the number in 'phi<number>.txt'.
    Return -1 if no number found to keep these files at the start.
    """
    match = re.search(r'phi(\d+)\.txt', value)
    if match:
        return int(match.group(1))
    else:
        return -1

# List and sort files by numerical index
phi_files = sorted(glob.glob(file_pattern), key=numerical_sort)

# Load coordinate arrays for x and y
coords_x = np.loadtxt(os.path.join(folder, 'coords_x.txt'))
coords_y = np.loadtxt(os.path.join(folder, 'coords_y.txt'))

# Subset of frames to animate
step_count_start = 0 
#step_count_end = 1500
frames = phi_files[step_count_start:]

# Create figure and axis for plotting
fig, ax = plt.subplots(figsize=(1, 1))

# Load and reshape first frame data for initialization
phi_1d = np.loadtxt(frames[0])
phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))

# Define custom palette colors as hex codes (blue, green, red)
palette = [ '#000fff', '#90ff70', '#ee0000']


'''
def get_boundaries(data):
    """
    Calculate boundaries for discretizing colormap based on data min and max.
    Adjust ratios as needed to restrict red to peak and blue to tail.
    """
    vmin, vmax = data.min(), data.max()
    return [vmin, vmin + 0.03*(vmax - vmin), vmin + 0.7*(vmax - vmin), vmax]
'''

color_steps = 25 
def get_boundaries(data, steps=10):
    vmin, vmax = data.min(), data.max()
    blue_end = vmin + 0.1 * (vmax - vmin)
    red_start = vmin + 0.95 * (vmax - vmin)
    # The boundaries for the 10 steps
    boundaries = np.linspace(blue_end, red_start, steps)
    # Complete set: minimum, blue_end, steps in-between, red_start, maximum
    boundaries = np.concatenate(([vmin, blue_end], boundaries, [red_start, vmax]))
    return boundaries

# Create segmented colormap from palette
segmented_cmap = LinearSegmentedColormap.from_list('red_green_blue', palette, N=color_steps)

# Initial boundaries and normalization
'''
boundaries = get_boundaries(phi)
norm = BoundaryNorm(boundaries, ncolors=segmented_cmap.N, clip=True)
'''
boundaries = get_boundaries(phi, color_steps)
num_color_intervals = len(boundaries) -1
norm = BoundaryNorm(boundaries, ncolors=num_color_intervals, clip=True)

# Plot initial heatmap with custom colormap and norm
im = ax.imshow(phi,
               extent=[coords_x[0], coords_x[-1], coords_y[0], coords_y[-1]],
               origin='lower',
               cmap=segmented_cmap,
               norm=norm,
               aspect='auto')  # Rectangular; use 'equal' if square pixels required

# Add colorbar with label
cbar = fig.colorbar(im, ax=ax, label='phi')

# Axis labels and aspect ratio
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal')

# Title for the initial frame
#title_text = ax.set_title(f'Phi heatmap from {os.path.basename(frames[0])}')
title_text = ax.text(0.5, 1.02, f'Phi heatmap at time = {time_stamps.iloc[step_count_start]}', transform=ax.transAxes, ha='center')

# Draw fixed horizontal line at y=0.3
Ys = 0.3
hline = ax.axhline(y=Ys, color='black', linewidth=2)

# Positions, colors, and widths for vertical lines (partial height)
v_positions = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]
colors_lines = ['gray', 'black'] * (len(v_positions)//2 + 1)
widths = [1, 2] * (len(v_positions)//2 + 1)

# Draw vertical lines from y=0 to y=Ys
v_lines = []
for x, c, w in zip(v_positions, colors_lines, widths):
    line = ax.axvline(x=x, ymin=0.0, ymax=Ys, color=c, linewidth=w)
    v_lines.append(line)

def update(frame_idx):
    """
    Update function for each animation frame:
    - Load and reshape frame data
    - Recalculate boundaries and normalization
    - Update image data and normalization
    - Update title
    """
    filename = frames[frame_idx]
    phi_1d = np.loadtxt(filename)
    phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))

    # Recalculate boundaries and norm
    boundaries = get_boundaries(phi,color_steps)
    num_color_intervals = len(boundaries) -1
    norm = BoundaryNorm(boundaries, ncolors=num_color_intervals, clip=True)
    im.set_norm(norm)

    # Update image data (transpose to match orientation)
    im.set_data(phi)

    # Update title text
    #title_text.set_text(f'Phi heatmap from {os.path.basename(filename)}')
    title_text.set_text(f'Phi heatmap at time = {time_stamps.iloc[frame_idx]}')


    return [im, title_text]

# Create animation using FuncAnimation
ani = FuncAnimation(fig, update, frames=len(frames), interval=100, blit=False, repeat=False)

# Show the plot and animation
plt.show()

