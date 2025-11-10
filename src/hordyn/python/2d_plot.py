import numpy as np
import matplotlib.pyplot as plt


edges_x= np.loadtxt('output/f0/coords_x.txt')
edges_y = np.loadtxt('output/f0/coords_y.txt')
phi_1d = np.loadtxt('output/f0/phi0.txt')

phi = phi_1d.reshape((len(edges_y)-1, len(edges_x)-1))

X, Y = np.meshgrid(edges_x, edges_y)

plt.figure(figsize=(8, 1.0))

# Set colormap to blue-white-red
pcm = plt.pcolormesh(X, Y, phi, shading='auto', cmap='bwr')
plt.colorbar(pcm, label='phi')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Phi[x,y] heatmap')

ax = plt.gca()
ax.set_aspect('auto')  # Or set aspect ratio if needed

plt.show()
