import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob
import os
import re
import pandas as pd

base_folder = 'output/hordyn-2016/'
folders = ['f0']
num_cycles = 1

def numerical_sort(value):
    match = re.search(r'phi(\d+)\.txt', value)
    return int(match.group(1)) if match else -1

coords_x = np.loadtxt(os.path.join(base_folder, 'f0', 'coords_x.txt'))
coords_y = np.loadtxt(os.path.join(base_folder, 'f0', 'coords_y.txt'))

fig, axs = plt.subplots(nrows=len(folders), ncols=1, figsize=(8, 20), constrained_layout=True)
if len(folders) == 1:
    axs = [axs]

im_list = []
phi_files_list = []
time_stamps_list = []

# Preload all phi files and time stamps for each folder
for idx, folder_name in enumerate(folders):
    folder = os.path.join(base_folder, folder_name)
    phi_files = sorted(glob.glob(os.path.join(folder, 'phi*.txt')), key=numerical_sort)
    phi_files_list.append(phi_files)
    df = pd.read_csv(os.path.join(folder, 'summary.csv'))
    time_stamps_list.append(df['time'])

    # Load first frame to initialize the image plot
    phi_1d = np.loadtxt(phi_files[0])
    phi = phi_1d.reshape((len(coords_y) - 1, len(coords_x) - 1))

    im = axs[idx].imshow(phi,
                        extent=[coords_x[0], coords_x[-1], coords_y[0], coords_y[-1]],
                        origin='lower',
                        aspect='equal')  # Default colormap and normalization
    axs[idx].set_ylabel(folder_name)
    axs[idx].set_aspect('equal')
    axs[idx].axhline(y=0.5, color='black', linewidth=2)

    # Additional vertical lines if required
    v_positions = np.arange(0.5, num_cycles, 0.5)
    colors_lines = ['gray', 'black'] * (len(v_positions) // 2 + 1)
    widths = [1, 2] * (len(v_positions) // 2 + 1)
    for x, c, w in zip(v_positions, colors_lines, widths):
        axs[idx].axvline(x=x, ymin=0, ymax=0.3, color=c, linewidth=w)

    if idx < 9:
        axs[idx].set_xticklabels([])
    else:
        axs[idx].set_xlabel('x')

    im_list.append(im)

# Create a single colorbar for the first image and link it to all axes
cbar = fig.colorbar(im_list[0], ax=axs, orientation='vertical')
cbar.set_label('Phi value')

title_texts = []
for idx in range(len(folders)):
    # Create text inside the axes, near the top center for frame info
    t = axs[idx].text(0.5, 1.05, '', transform=axs[idx].transAxes, ha='center')
    title_texts.append(t)

def update(frame_idx):
    artists = []
    for idx in range(len(folders)):
        phi_1d = np.loadtxt(phi_files_list[idx][frame_idx])
        phi = phi_1d.reshape((len(coords_y) - 1, len(coords_x) - 1))

        im_list[idx].set_data(phi)
        # No normalization update needed, keep fixed colormap scaling

        title_texts[idx].set_text(f'{folders[idx]} at time = {time_stamps_list[idx].iloc[frame_idx]}')

        artists.append(im_list[idx])
        artists.append(title_texts[idx])
    return artists

ani = FuncAnimation(fig, update, frames=len(phi_files_list[0]), interval=100, blit=False, repeat=False)

plt.show()

