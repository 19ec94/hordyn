import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

Dc = 1.0
x_s = 0.5
g_l = 0.5  # velocity for x < x_s
g_r = 1.0  # velocity for x >= x_s

def phi_0(x):
    return np.where(x <= x_s, 1.0, 0.5)

def phi(t, x):
    # Calculate characteristic foot points going backward in time
    x_char_l = x - g_l * t
    x_char_r = x - g_r * t

    # Count number of domain crossings (how many times flux doubled)
    n_cross_l = np.floor(x_char_l / Dc) - np.floor(x / Dc)
    n_cross_r = np.floor(x_char_r / Dc) - np.floor(x / Dc)

    # Map foot points back to domain [0, Dc)
    x_mod_l = np.mod(x_char_l, Dc)
    x_mod_r = np.mod(x_char_r, Dc)

    # Select which velocity regime applies based on foot point
    cond = x_mod_l < x_s

    # Flux doubling factor from number of crossings (negative sign because characteristics traced backward)
    factor_l = 2 ** (-n_cross_l)
    factor_r = 2 ** (-n_cross_r)

    return np.where(cond, phi_0(x_mod_l) * factor_l, phi_0(x_mod_r) * factor_r)

x = np.linspace(0, Dc, 1000)

fig, ax = plt.subplots()
line, = ax.plot(x, phi(0, x))
ax.set_ylim(0, 5)
ax.set_xlabel('x')
ax.set_ylabel(r'$\phi(t,x)$')
ax.set_title('Analytical solution with flux doubling animation')

def update(frame):
    t = frame * 0.01  # time step
    y = phi(t, x)
    line.set_ydata(y)
    ax.set_title(f'Time={t:.2f}')
    return line,

ani = FuncAnimation(fig, update, frames=100, interval=100, blit=False, repeat=False)

plt.show()

