import numpy as np
import matplotlib.pyplot as plt

Dc = 1.0
x = np.linspace(0, 1, 500)  # 500 points between 0 and 1
y = 2 ** (-x / Dc)

plt.plot(x, y)
plt.xlabel('x')
plt.ylabel(r'$2^{-x/D_c}$')
plt.title('Plot of $2^{-x/D_c}$ with $D_c=1.0$')
plt.grid(True)
plt.show()

