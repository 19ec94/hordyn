import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import re

# Directory containing CSV files
data_dir = './2dplain_results'
Ys = 0.3;
# Get list of all files in directory matching pattern 'step_*.csv'
all_files = [f for f in os.listdir(data_dir) if re.match(r'step_\d+\.csv', f)]

# Extract step numbers and sort files by these numbers
def extract_step_number(filename):
    match = re.search(r'step_(\d+)\.csv', filename)
    return int(match.group(1)) if match else -1

sorted_files = sorted(all_files, key=extract_step_number)

# Full paths
file_paths = [os.path.join(data_dir, f) for f in sorted_files]

num_frames = len(file_paths)  # Total number of frames from your files

# Select about 50 frame indices evenly spaced throughout the animation
num_save_frames = 50
save_frame_indices = set(np.linspace(0, num_frames - 1, num_save_frames, dtype=int))

# Read the first CSV file to get coordinates etc.
initial_data = pd.read_csv(file_paths[0])
x_unique = np.sort(initial_data['x'].unique())
y_unique = np.sort(initial_data['y'].unique())
Nx = x_unique.size
Ny = y_unique.size

# Plot setup
fig, ax = plt.subplots(figsize=(6,5))
phi = initial_data['phi'].values.reshape((Nx, Ny))
im = ax.imshow(phi.T, extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
               origin='lower', cmap='bwr', vmin=phi.min(), vmax=phi.max())
title_text = ax.set_title('Full biological model, t=0', fontsize=15)
ax.set_xlabel('$a$', fontsize=20)
ax.set_ylabel('$\gamma$', fontsize=20)
hline = ax.axhline(y=Ys, color='black', linewidth=2)
v1 = ax.axvline(x=0.5, ymin=0.0, ymax=Ys,  color='gray', linewidth=1)
v2 = ax.axvline(x=1.0, ymin=0.0, ymax=Ys,  color='black', linewidth=2)
v3 = ax.axvline(x=1.5, ymin=0.0, ymax=Ys,  color='gray', linewidth=1)
v4 = ax.axvline(x=2.0, ymin=0.0, ymax=Ys,  color='black', linewidth=2)
v5 = ax.axvline(x=2.5, ymin=0.0, ymax=Ys,  color='gray', linewidth=1)
v6 = ax.axvline(x=3.0, ymin=0.0, ymax=Ys,  color='black', linewidth=2)
v7 = ax.axvline(x=3.5, ymin=0.0, ymax=Ys,  color='gray', linewidth=1)
v8 = ax.axvline(x=4.0, ymin=0.0, ymax=Ys,  color='black', linewidth=2)
v9 = ax.axvline(x=4.5, ymin=0.0, ymax=Ys,  color='gray', linewidth=1)
v10 = ax.axvline(x=5.0, ymin=0.0, ymax=Ys,  color='black', linewidth=2)
v11 = ax.axvline(x=5.5, ymin=0.0, ymax=Ys,  color='gray', linewidth=1)
v12 = ax.axvline(x=6.0, ymin=0.0, ymax=Ys,  color='black', linewidth=2)
v11 = ax.axvline(x=6.5, ymin=0.0, ymax=Ys,  color='gray', linewidth=1)
v12 = ax.axvline(x=7.0, ymin=0.0, ymax=Ys,  color='black', linewidth=2)
v11 = ax.axvline(x=7.5, ymin=0.0, ymax=Ys,  color='gray', linewidth=1)

# You may want to keep references to all vertical lines here
v_lines = [v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11]  # add all lines you declared

saved_frame_count = 0  # counter for saved frames

def update(frame):
    global saved_frame_count
    filename = file_paths[frame]
    data = pd.read_csv(filename)
    phi = data['phi'].values.reshape((Nx, Ny))
    im.set_data(phi.T)

    # Extract step number to show in title
    step_num = extract_step_number(os.path.basename(filename))

    # Your text annotations - consider creating them once outside update or clear them each frame if needed
    ax.text(0.1, 0.1, r'$\Omega_1^1$', color='white', fontsize=12, weight='bold')
    ax.text(0.6, 0.1, r'$\Omega_2^1$', color='white', fontsize=12, weight='bold')
    ax.text(1.1, 0.1, r'$\Omega_1^2$', color='white', fontsize=12, weight='bold')
    ax.text(1.6, 0.1, r'$\Omega_2^2$', color='white', fontsize=12, weight='bold')
    ax.text(2.1, 0.1, r'$\Omega_1^3$', color='white', fontsize=12, weight='bold')
    ax.text(2.6, 0.1, r'$\Omega_2^3$', color='white', fontsize=12, weight='bold')
    ax.text(3.1, 0.1, r'$\Omega_1^4$', color='white', fontsize=12, weight='bold')
    ax.text(3.6, 0.1, r'$\Omega_2^4$', color='white', fontsize=12, weight='bold')
    ax.text(4.1, 0.1, r'$\Omega_1^5$', color='white', fontsize=12, weight='bold')
    ax.text(4.6, 0.1, r'$\Omega_2^5$', color='white', fontsize=12, weight='bold')
    ax.text(5.1, 0.1, r'$\Omega_1^6$', color='white', fontsize=12, weight='bold')
    ax.text(5.6, 0.1, r'$\Omega_2^6$', color='white', fontsize=12, weight='bold')
    ax.text(6.1, 0.1, r'$\Omega_1^7$', color='white', fontsize=12, weight='bold')
    ax.text(6.6, 0.1, r'$\Omega_2^7$', color='white', fontsize=12, weight='bold')
    ax.text(7.1, 0.1, r'$\Omega_1^8$', color='white', fontsize=12, weight='bold')
    ax.text(7.6, 0.1, r'$\Omega_2^8$', color='white', fontsize=12, weight='bold')
    ax.text(4.0, 0.8, r'$\Omega_3$', color='white', fontsize=12, weight='bold')

    title_text.set_text(f'Transport of $\phi(a, \gamma)$ at step={step_num}')

    # Save only if it's one of the selected frames
    #if frame in save_frame_indices:
    #    saved_frame_count += 1
    #    filename_save = os.path.join(data_dir, f"frame_{saved_frame_count}.png")
    #    fig.savefig(filename_save, bbox_inches='tight')
    #    print(f"Saved {filename_save}")

    # Return all artists that get updated
    return [im, hline] + v_lines + [title_text]

ani = FuncAnimation(fig, update, frames=num_frames, interval=200, blit=False, repeat=False)

plt.show()
