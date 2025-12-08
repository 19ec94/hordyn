import numpy as np
import matplotlib.pyplot as plt

Dc = 1.0
x_s = 0.5
g_l = 0.5  # velocity for x < 0.5
g_r = 1.0  # velocity for x >= 0.5

def phi_0(x):
    return np.where(x < x_s, 1.0, 0.5)

def phi(t, x):
    # Compute foot points of characteristics
    x_char_l = x - g_l * t
    x_char_r = x - g_r * t

    # Count how many times characteristics crossed boundary x=1
    n_cross_l = np.floor(x_char_l / Dc) - np.floor(x / Dc)
    n_cross_r = np.floor(x_char_r / Dc) - np.floor(x / Dc)

    # Map foot points back into [0, Dc)
    x_mod_l = np.mod(x_char_l, Dc)
    x_mod_r = np.mod(x_char_r, Dc)

    condition = x_mod_l < x_s

    # Flux doubling factor: 2^number_crossings
    factor_l = 2 ** (-n_cross_l)  # negative because of backwards char
    factor_r = 2 ** (-n_cross_r)

    # Compute solution from initial condition and flux doubling
    val_l = phi_0(x_mod_l) * factor_l
    val_r = phi_0(x_mod_r) * factor_r

    return np.where(condition, val_l, val_r)


x = np.linspace(0, Dc, 1000)
times = [0.0, 0.3, 0.6, 0.9, 1.2]

plt.figure(figsize=(8,5))
for t in times:
    plt.plot(x, phi(t, x), label=f't = {t:.1f}')

plt.xlabel('x')
plt.ylabel(r'$\phi(t,x)$')
plt.title('Analytical solution with flux doubling at $x=1$ and periodic BC')
plt.legend()
plt.grid(True)
plt.show()

