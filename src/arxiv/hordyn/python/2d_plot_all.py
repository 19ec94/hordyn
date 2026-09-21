import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import glob
import os
import re

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

fig, ax = plt.subplots(figsize=(8,1))
pcm = None # to hold the plot object

# Loop through the desired number of files and plot
for i, filepath in enumerate(phi_files[step_count_start:step_count_end]):
    phi_1d = np.loadtxt(filepath)
    phi = phi_1d.reshape((len(coords_y)-1, len(coords_x)-1))

    if pcm is None:
        #pcm = ax.pcolormesh(X, Y, phi, shading='auto', cmap='bwr')
        pcm = ax.pcolormesh(X, Y, phi, shading='auto', cmap='bwr', norm=colors.SymLogNorm( linthresh=0.03, linscale=0.03,vmin=phi.min(), vmax=phi.max(), base=10))
        fig.colorbar(pcm, ax=ax, label='phi')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    else:
        pcm.set_array(phi.ravel()) #update data
        pcm.set_clim(np.min(phi), np.max(phi))  #update color limits if necessary
    ax.set_title(f'Phi heatmap from {os.path.basename(filepath)}')
    plt.pause(0.1)
plt.show()
