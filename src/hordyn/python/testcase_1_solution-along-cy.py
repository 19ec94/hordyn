import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re
from matplotlib.animation import FuncAnimation

folder = 'output/1d_testcase_1/gl2.0/f0'
file_pattern = os.path.join(folder, 'phi*.txt')
centers_x = np.loadtxt(os.path.join(folder, 'centers_x.txt'))
centers_y = np.loadtxt(os.path.join(folder, 'centers_y.txt'))

# Read the timestamps
import pandas as pd
df = pd.read_csv(os.path.join(folder, 'summary.csv'))
time_stamps = df['time'] 

num_cols = len(centers_y)
ly = 1.0
dy = ly / num_cols

y_s = 0.15
index = int(y_s / dy)

def numerical_sort(value):
    match = re.search(r'phi(\d+)\.txt', value)
    return int(match.group(1)) if match else -1

phi_files = sorted(glob.glob(file_pattern), key=numerical_sort)
step_count_start = 0 
step_count_end = 10
files_to_plot = phi_files[step_count_start:]

fig, ax = plt.subplots(figsize=(1,1))
line, = ax.plot([], [], label=f'Row at y = {y_s}')  # empty line initially

ax.set_xlabel('x')
ax.set_ylabel('phi')
#ax.set_title(f'Row data at y = {y_s} (index {index})')
ax.set_title(f'Phi heatmap at time = {time_stamps.iloc[step_count_start]}')
ax.set_xlim(centers_x.min(), centers_x.max())

def init():
    line.set_data([], [])
    return line,

def update(frame_idx):
    phi_file = files_to_plot[frame_idx]
    with open(phi_file, 'r') as file:
        for i, row in enumerate(file):
            if i == index:
                row_data = np.fromstring(row, sep=' ')
                break
    line.set_data(centers_x, row_data)
    #ax.set_title(f'{os.path.basename(phi_file)} at y = {y_s}')
    ax.set_title(f'Phi heatmap at time = {time_stamps.iloc[frame_idx]}')
    ymin, ymax = row_data.min(), row_data.max()
    ax.set_ylim(ymin, ymax)
    return line,

ani = FuncAnimation(fig, update, frames=len(files_to_plot), init_func=init, blit=False, interval=500, repeat=False)
plt.legend()
plt.show()

