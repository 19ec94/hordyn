import pandas as pd
import matplotlib.pyplot as plt

dfs = []
for i in range(1):
    filename = f'output/hordyn-2016/f{i}/summary.csv'
    df = pd.read_csv(filename)
    dfs.append(df)

fig, ax1 = plt.subplots(figsize=(8, 6), constrained_layout=True)

for i, df in enumerate(dfs):
    ax1.plot(df['time'], df['mass_mitosis']/df['mass'], label=f"mass f{i}")

ax1.set_xlabel('Time')
ax1.set_ylabel('mass')
ax1.grid(True)

# Create a second y-axis sharing the same x-axis
ax2 = ax1.twinx()

for i, df in enumerate(dfs):
    ax2.plot(df['time'], df['mass_mitosis']/df['mass'], linestyle='--', label=f"MI f{i}")

ax2.set_ylabel('MI')
#ax2.set_ylim(0.0016, 0.0056)

# Combine legends from both axes
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()


#fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 6), constrained_layout=True)
'''
variable = "global_maturity"
for i, df in enumerate(dfs):
    axs[0,0].plot(df['time'], df[variable], label=f"f{i}")

axs[0,0].set_xlabel('Time')
axs[0,0].set_ylabel(variable)
axs[0,0].set_title(f'Time vs {variable}')
axs[0,0].grid(True)
axs[0,0].axhline(y=5.0, linestyle='--', linewidth=2)
axs[0,0].legend()

variable = "MI"
for i, df in enumerate(dfs):
    axs[0,1].plot(df['time'], df[variable], label=f"f{i}")

axs[0,1].set_xlabel('Time')
axs[0,1].set_ylabel(variable)
axs[0,1].set_title(f'Time vs {variable}')
axs[0,1].grid(True)
axs[0,1].legend()

variable = "phi_max"
for i, df in enumerate(dfs):
    axs[1,0].plot(df['time'], df[variable], label=f"f{i}")

axs[1,0].set_xlabel('Time')
axs[1,0].set_ylabel(variable)
axs[1,0].set_title(f'Time vs {variable}')
axs[1,0].grid(True)
axs[1,0].legend()

variable = "mass"
for i, df in enumerate(dfs):
    axs[1,1].plot(df['time'], df[variable], label=f"f{i}")

axs[1,1].set_xlabel('Time')
axs[1,1].set_ylabel(variable)
axs[1,1].set_title(f'Time vs {variable}')
axs[1,1].grid(True)
axs[1,1].legend()
'''


'''
variable = "global_fsh"
for i, df in enumerate(dfs):
    axs[1,0].plot(df['time'], df[variable], label=f"f{i}")

axs[1,0].set_xlabel('Time')
axs[1,0].set_ylabel(variable)
axs[1,0].set_title(f'Time vs {variable}')
axs[1,0].grid(True)
axs[1,0].legend()

variable = "local_fsh"
for i, df in enumerate(dfs):
    axs[1,1].plot(df['time'], df[variable], label=f"f{i}")

axs[1,1].set_xlabel('Time')
axs[1,1].set_ylabel(variable)
axs[1,1].set_title(f'Time vs {variable}')
axs[1,1].grid(True)
axs[1,1].legend()
'''

plt.show()

