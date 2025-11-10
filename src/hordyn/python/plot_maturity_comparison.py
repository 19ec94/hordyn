import pandas as pd
import matplotlib.pyplot as plt
import os

data_path_1 = 'output/fullmodel/2012_tau1.0/f0/summary.csv'
df_1 = pd.read_csv(data_path_1)
output_dir_1 = os.path.dirname(data_path_1)

data_path_2 = 'output/fullmodel/2012_tau1.2/f0/summary.csv'
df_2 = pd.read_csv(data_path_2)
output_dir_2 = os.path.dirname(data_path_2)

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

if show_plot_pairs:
    # Define pairs of indices in plot_vars for side-by-side plotting
    pairs = [
            #(0, 1),  # phi_max & mass
            (2, 3),  # local_maturity & global_maturity
            #(4, 5)   # local_fsh & global_fsh
            ]

    for i, (idx1, idx2) in enumerate(pairs):
        fig, axs = plt.subplots(1, 2, figsize=(16, 4))

        # Plot first variable on left subplot
        v1 = plot_vars[idx1]
        ax = axs[0]
        if v1['color']:
            ax.plot(df_1['time'], df_1[v1['column']], label="tau=1.0")
            ax.plot(df_2['time'], df_2[v1['column']], label="tau=1.2")
        else:
            ax.plot(df_1['time'], df_1[v1['column']], label="tau=1.0")
            ax.plot(df_2['time'], df_2[v1['column']], label="tau=1.2")
        ax.set_xlabel('Time')
        ax.set_ylabel(v1['ylabel'])
        ax.set_title(v1['title'])
        ax.grid(True)
        ax.legend()

        # Plot second variable on right subplot
        v2 = plot_vars[idx2]
        ax = axs[1]
        if v2['color']:
            ax.plot(df_1['time'], df_1[v2['column']], label="tau=1.0")
            ax.plot(df_2['time'], df_2[v2['column']], label="tau=1.2")
        else:
            ax.plot(df_1['time'], df_1[v2['column']], label="tau=1.0")
            ax.plot(df_2['time'], df_2[v2['column']], label="tau=1.2")
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

