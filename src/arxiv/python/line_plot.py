import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV data
data = pd.read_csv('2dplain_results/step_1200.csv')

# Target y value and tolerance for filtering
y_target = 0.2
tolerance = 1e-2

# Filter data close to y=0.15
subset = data[(data['y'] >= y_target - tolerance) & (data['y'] <= y_target + tolerance)]

# Sort by x for an ordered plot
subset = subset.sort_values('x')

# Plot
plt.plot(subset['x'], subset['phi'])
plt.xlabel('x')
plt.ylabel(r'$\phi(x, 0.25)$')
plt.title('Plot of φ along y=0.25')
plt.grid(True)
plt.show()

