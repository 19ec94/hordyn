import numpy as np
import matplotlib.pyplot as plt

def read_grid_data(filename):
    with open(filename, 'r') as f:
        sizes = list(map(int, f.readline().split()))
        Nx1, Ny1, Nx, Ny, phiNx, phiNy = sizes

        xCoords = np.array(list(map(float, f.readline().split())))
        yCoords = np.array(list(map(float, f.readline().split())))
        xCenters = np.array(list(map(float, f.readline().split())))
        yCenters = np.array(list(map(float, f.readline().split())))

        phi_data = []
        for _ in range(phiNx):
            row = list(map(float, f.readline().split()))
            phi_data.append(row)
        phi = np.array(phi_data)

    return xCoords, yCoords, xCenters, yCenters, phi

def plot_grid_and_phi(xCoords, yCoords, xCenters, yCenters, phi):
    # Create meshgrid for cell centers where phi is defined
    Xc, Yc = np.meshgrid(xCenters, yCenters, indexing='ij')

    # Plot phi as a heatmap over cell centers
    plt.figure(figsize=(8, 5))
    plt.pcolormesh(xCoords, yCoords, phi.T, shading='auto', cmap='viridis')
    plt.colorbar(label='Phi value')

    # Overlay grid lines from edges
    for x in xCoords:
        plt.axvline(x=x, color='k', linewidth=0.5, linestyle='--')
    for y in yCoords:
        plt.axhline(y=y, color='k', linewidth=0.5, linestyle='--')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Phi values with grid overlay')
    plt.tight_layout()
    plt.show()


# Usage:
xCoords, yCoords, xCenters, yCenters, phi = read_grid_data('../../InitialData.txt')
plot_grid_and_phi(xCoords, yCoords, xCenters, yCenters, phi)

