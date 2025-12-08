import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import glob
import re

folder = 'output/hordyn-2016/f0'
centers_x = np.loadtxt(os.path.join(folder, 'centers_x.txt'))

index = 0 

def read_line(file_path, line_num):
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            if i == line_num:
                return np.fromstring(line, sep=' ')
    raise ValueError(f"Line {line_num} not found in file.")

def numerical_sort(value):
    match = re.search(r'phi(\d+)\.txt', value)
    return int(match.group(1)) if match else -1

phi_files = sorted(glob.glob(os.path.join(folder, 'phi*.txt')), key=numerical_sort)

# For debuggling
y1 = read_line(phi_files[0], index)
y2 = read_line(phi_files[-1], index)
plt.plot(centers_x, y1)
plt.plot(centers_x, y2)
plt.show()

'''
fig, ax = plt.subplots()
line, = ax.plot(centers_x, np.zeros_like(centers_x))
ax.set_ylim(0, 8)

def init():
    line.set_ydata(np.zeros_like(centers_x))
    return line,

def animate(frame):
    row_data = read_line(phi_files[frame], index)
    #print(row_data) # correctly prints the data
    line.set_ydata(row_data)
    ax.relim()
    ax.autoscale_view()
    t = frame * 0.002
    ax.set_title(f"Time: {t}") # correctly prints the frame number
    return line,

ani = FuncAnimation(fig, animate, frames=len(phi_files), init_func=init,
                    interval=100, blit=False, repeat=False)
plt.show()
'''
