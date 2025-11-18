import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from matplotlib.animation import FuncAnimation
import glob
import os
import re
import pandas as pd

base_folder = 'output/fullmodel/parallel/10follicles'
#folders = [f'f{i}' for i in range(3)]
folders = ['f0', 'f1', 'f2', 'f3', 'f4', 'f9']

def numerical_sort(value):
    match = re.search(r'phi(\d+)\.txt', value)
    return int(match.group(1)) if match else -1

coords_x = np.loadtxt(os.path.join(base_folder, 'f0', 'coords_x.txt'))
coords_y = np.loadtxt(os.path.join(base_folder, 'f0', 'coords_y.txt'))

palette = ['#000fff', '#90ff70', '#ee0000']
color_steps = 25
segmented_cmap = LinearSegmentedColormap.from_list('red_green_blue', palette, N=color_steps)

def get_boundaries(data, steps=color_steps):
    vmin, vmax = data.min(), data.max()
    blue_end = vmin + 0.1 * (vmax - vmin)
    red_start = vmin + 0.95 * (vmax - vmin)
    boundaries = np.linspace(blue_end, red_start, steps)
    boundaries = np.concatenate(([vmin, blue_end], boundaries, [red_start, vmax]))
    return boundaries

fig, axs = plt.subplots(nrows=len(folders), ncols=1, figsize=(8, 20), constrained_layout=True)

im_list = []
phi_files_list = []
time_stamps_list = []

# Initialize all images once
for idx, folder_name in enumerate(folders):
    folder = os.path.join(base_folder, folder_name)
    phi_files = sorted(glob.glob(os.path.join(folder, 'phi*.txt')), key=numerical_sort)
    phi_files_list.append(phi_files)
    df = pd.read_csv(os.path.join(folder, 'summary.csv'))
    time_stamps_list.append(df['time'])

    phi_1d = np.loadtxt(phi_files[0])
    phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))

    boundaries = get_boundaries(phi, color_steps)
    norm = BoundaryNorm(boundaries, ncolors=len(boundaries)-1, clip=True)

    im = axs[idx].imshow(phi,
                        extent=[coords_x[0], coords_x[-1], coords_y[0], coords_y[-1]],
                        origin='lower',
                        cmap=segmented_cmap,
                        norm=norm,
                        aspect='equal')
    axs[idx].set_ylabel(folder_name)
    axs[idx].set_aspect('equal')
    axs[idx].axhline(y=0.3, color='black', linewidth=2)

    v_positions = np.arange(0.5, 4.0, 0.5)
    colors_lines = ['gray', 'black'] * (len(v_positions)//2 + 1)
    widths = [1, 2] * (len(v_positions)//2 + 1)
    for x, c, w in zip(v_positions, colors_lines, widths):
        axs[idx].axvline(x=x, ymin=0, ymax=0.3, color=c, linewidth=w)

    if idx < 9:
        axs[idx].set_xticklabels([])
    else:
        axs[idx].set_xlabel('x')

    im_list.append(im)

title_texts = []
for idx in range(len(folders)):
    # Create text inside the axes, near top center
    t = axs[idx].text(0.5, 1.05, '', transform=axs[idx].transAxes, ha='center')
    title_texts.append(t)

def update(frame_idx):
    artists = []
    for idx in range(len(folders)):
        phi_1d = np.loadtxt(phi_files_list[idx][frame_idx])
        phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))

        boundaries = get_boundaries(phi, color_steps)
        norm = BoundaryNorm(boundaries, ncolors=len(boundaries)-1, clip=True)

        im_list[idx].set_data(phi)
        im_list[idx].set_norm(norm)

        # Update title text instead of set_title
        title_texts[idx].set_text(f'{folders[idx]} at time = {time_stamps_list[idx].iloc[frame_idx]}')

        artists.append(im_list[idx])
        artists.append(title_texts[idx])
    return artists

ani = FuncAnimation(fig, update, frames=len(phi_files_list[0]), interval=100, blit=False, repeat=False)

plt.show()

