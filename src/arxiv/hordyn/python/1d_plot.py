import numpy as np
import matplotlib.pyplot as plt
import os

folder = 'output/1d_testcase_1/gl2.0/f0'
file_path = os.path.join(folder, 'phi500.txt')

centers_x = np.loadtxt(os.path.join(folder, 'centers_x.txt'))

num_cols = 200
ly = 1.0
dy = ly / num_cols

# y position of interest
y_s = 0.15

# Compute row index corresponding to y_s
index = int(y_s / dy)
print(f"Reading line index: {index}")

# Function to read specific line (0-indexed) from a file
def read_line(file_path, line_num):
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            if i == line_num:
                # Convert line string to numpy array of floats
                return np.fromstring(line, sep=' ')
    raise ValueError(f"Line {line_num} not found in file.")

# Read only the specific line
row_data = read_line(file_path, index)

# Plot the row data versus centers_x
plt.plot(centers_x, row_data, label=f'Row at y index {index}')
plt.xlabel('x')
plt.ylabel('phi')
plt.title(f'Row data at y = {y_s} (index {index})')
plt.legend()
plt.show()

