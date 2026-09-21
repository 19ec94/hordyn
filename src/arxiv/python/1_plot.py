import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re

# Directory containing CSV files
data_dir = './2dplain_results'
Ys = 0.3

# Function to extract step number from filename
def extract_step_number(filename):
    match = re.search(r'step_(\d+)\.csv', filename)
    return int(match.group(1)) if match else -1

# List and sort files by step number
all_files = [f for f in os.listdir(data_dir) if re.match(r'step_\d+\.csv', f)]
sorted_files = sorted(all_files, key=extract_step_number)
file_paths = [os.path.join(data_dir, f) for f in sorted_files]

# Choose a single file to plot (e.g., first file)
file_to_plot = file_paths[1550]

# Load data from the chosen file
data = pd.read_csv(file_to_plot)
x_unique = np.sort(data['x'].unique())
y_unique = np.sort(data['y'].unique())
Nx = x_unique.size
Ny = y_unique.size

# Reshape phi data for plotting
phi = data['phi'].values.reshape((Nx, Ny))

# Plot setup
fig, ax = plt.subplots(figsize=(6,5))
im = ax.imshow(phi.T, extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
               origin='lower', cmap='bwr', vmin=phi.min(), vmax=phi.max())

# Add title and labels
step_num = extract_step_number(os.path.basename(file_to_plot))
ax.set_title(f'Transport of $\phi(a, \gamma)$ at step={step_num}', fontsize=15)
ax.set_xlabel('$a$', fontsize=20)
ax.set_ylabel('$\\gamma$', fontsize=20)

# Add lines
ax.axhline(y=Ys, color='black', linewidth=2)
for x_pos, color, width in zip([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5],
                               ['gray', 'black', 'gray', 'black', 'gray', 'black', 'gray', 'black', 'gray', 'black', 'gray', 'black', 'gray', 'black', 'gray'],
                               [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1]):
    ax.axvline(x=x_pos, ymin=0.0, ymax=Ys, color=color, linewidth=width)

# Add text annotations once
texts = [
    (0.1, 0.1, r'$\Omega_1^1$'), (0.6, 0.1, r'$\Omega_2^1$'), (1.1, 0.1, r'$\Omega_1^2$'),
    (1.6, 0.1, r'$\Omega_2^2$'), (2.1, 0.1, r'$\Omega_1^3$'), (2.6, 0.1, r'$\Omega_2^3$'),
    (3.1, 0.1, r'$\Omega_1^4$'), (3.6, 0.1, r'$\Omega_2^4$'), (4.1, 0.1, r'$\Omega_1^5$'),
    (4.6, 0.1, r'$\Omega_2^5$'), (5.1, 0.1, r'$\Omega_1^6$'), (5.6, 0.1, r'$\Omega_2^6$'),
    (6.1, 0.1, r'$\Omega_1^7$'), (6.6, 0.1, r'$\Omega_2^7$'), (7.1, 0.1, r'$\Omega_1^8$'),
    (7.6, 0.1, r'$\Omega_2^8$'), (4.0, 0.8, r'$\Omega_3$')
]
for x, y, text in texts:
    ax.text(x, y, text, color='white', fontsize=12, weight='bold')

plt.show()

