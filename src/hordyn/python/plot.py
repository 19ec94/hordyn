import pandas as pd
import matplotlib.pyplot as plt
import os

#data_path = 'output/1d_testcase_2/gr0.5/f0/summary.csv'
data_path = 'output/fullmodel/2012_tau1.0/f0/summary.csv'
output_dir = os.path.dirname(data_path)
df = pd.read_csv(data_path)
show_plot = True # or False to save
show_separate_plots = False 
show_plot_pairs = True

# Define the variables to plot with their attributes
plot_vars = [
        {'column': 'phi_max', 'ylabel': 'phi_max', 'title': 'Time vs phi_max', 'color': None},
        {'column': 'mass', 'ylabel': 'Mass', 'title': 'Time vs Mass', 'color': 'green'},
        {'column': 'local_maturity', 'ylabel': 'Local Maturity', 'title': 'Time vs Local Maturity', 'color': 'red'},
        {'column': 'global_maturity', 'ylabel': 'Global Maturity', 'title': 'Time vs Global Maturity', 'color': 'blue'},
        {'column': 'local_fsh', 'ylabel': 'Local FSH', 'title': 'Time vs Local FSH', 'color': 'orange'},
        {'column': 'global_fsh', 'ylabel': 'Global FSH', 'title': 'Time vs Global FSH', 'color': 'purple'}
        ]

# single plots
if show_separate_plots:
    for var in plot_vars:
        fig, ax = plt.subplots(figsize=(8, 4))
        if var['color']:
            ax.plot(df['time'], df[var['column']], label=var['column'], color=var['color'])
        else:
            ax.plot(df['time'], df[var['column']], label=var['column'])
        ax.set_xlabel('Time')
        ax.set_ylabel(var['ylabel'])
        ax.set_title(var['title'])
        ax.grid(True)
        ax.legend()
        if show_plot:
            plt.show()
        else:
            filename = f"{var['column']}.png"
            save_path = os.path.join(output_dir, filename)
            fig.savefig(save_path)
            plt.close(fig)

if show_plot_pairs:
    # Define pairs of indices in plot_vars for side-by-side plotting
    pairs = [
            (0, 1),  # phi_max & mass
            (2, 3),  # local_maturity & global_maturity
            (4, 5)   # local_fsh & global_fsh
            ]

    for i, (idx1, idx2) in enumerate(pairs):
        fig, axs = plt.subplots(1, 2, figsize=(16, 4))

        # Plot first variable on left subplot
        v1 = plot_vars[idx1]
        ax = axs[0]
        if v1['color']:
            ax.plot(df['time'], df[v1['column']], label=v1['column'], color=v1['color'])
        else:
            ax.plot(df['time'], df[v1['column']], label=v1['column'])
        ax.set_xlabel('Time')
        ax.set_ylabel(v1['ylabel'])
        ax.set_title(v1['title'])
        ax.grid(True)
        ax.legend()

        # Plot second variable on right subplot
        v2 = plot_vars[idx2]
        ax = axs[1]
        if v2['color']:
            ax.plot(df['time'], df[v2['column']], label=v2['column'], color=v2['color'])
        else:
            ax.plot(df['time'], df[v2['column']], label=v2['column'])
        ax.set_xlabel('Time')
        ax.set_ylabel(v2['ylabel'])
        ax.set_title(v2['title'])
        #ax.set_ylim(bottom=0, top=2.0)
        #ax.set_xlim(right=0.4)
        ax.grid(True)
        ax.legend()

        plt.tight_layout()

        if show_plot:
            plt.show()
        else:
            filename = f'side_by_side_{v1["column"]}_{v2["column"]}.png'
            save_path = os.path.join(output_dir, filename)
            fig.savefig(save_path)
            plt.close(fig)



