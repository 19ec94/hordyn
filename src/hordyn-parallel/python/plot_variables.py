import pandas as pd
import matplotlib.pyplot as plt

dfs = []
for i in range(10):
    filename = f'output/fullmodel/parallel/10follicles/f{i}/summary.csv'
    df = pd.read_csv(filename)
    dfs.append(df)

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(8, 6), constrained_layout=True)

variable = "global_maturity"
for i, df in enumerate(dfs):
    axs[0,0].plot(df['time'], df[variable], label=f"f{i}")

axs[0,0].set_xlabel('Time')
axs[0,0].set_ylabel(variable)
axs[0,0].set_title(f'Time vs {variable}')
axs[0,0].grid(True)
axs[0,0].axhline(y=5.0, linestyle='--', linewidth=2)
axs[0,0].legend()

variable = "local_maturity"
for i, df in enumerate(dfs):
    axs[0,1].plot(df['time'], df[variable], label=f"f{i}")

axs[0,1].set_xlabel('Time')
axs[0,1].set_ylabel(variable)
axs[0,1].set_title(f'Time vs {variable}')
axs[0,1].grid(True)
axs[0,1].legend()

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

variable = "phi_max"
for i, df in enumerate(dfs):
    axs[2,0].plot(df['time'], df[variable], label=f"f{i}")

axs[2,0].set_xlabel('Time')
axs[2,0].set_ylabel(variable)
axs[2,0].set_title(f'Time vs {variable}')
axs[2,0].grid(True)
axs[2,0].legend()

variable = "mass"
for i, df in enumerate(dfs):
    axs[2,1].plot(df['time'], df[variable], label=f"f{i}")

axs[2,1].set_xlabel('Time')
axs[2,1].set_ylabel(variable)
axs[2,1].set_title(f'Time vs {variable}')
axs[2,1].grid(True)
axs[2,1].legend()

plt.show()

