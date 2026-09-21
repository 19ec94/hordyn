import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import PowerNorm


# ===================
# Global style settings
# ===================
# To plot COW data enter COW or PCOS.
case = "COW"

if (case == "COW"):
	dname = ["cow-beta0-0/", "cow-beta1-1/"]
	files = [ "phi0.csv", "phi1.csv", "phi2.csv"]
	ofname = "cow419-wave1-phi.pdf"
	nrow = 3
	ncol = 2
elif (case == "PCOS"):
	dname = ["pcos-beta0-0/", "pcos-beta1-3/"]
	files = [ "phi0.csv", "phi1.csv", "phi2.csv"]
	ofname = "pcos_heatmap.pdf"
	nrow=3
	ncol=2
elif (case == "AYMARD"):
	dname = ["aymard-tmp/"]
	files = [ "phi0.csv", "phi1.csv", "phi2.csv", "phi3.csv"]
	ofname = "aymard_heatmap.pdf"
	nrow=4
	ncol=2
elif (case == "AYMARD-2FOL"):
	dname = ["aymard-2fol-f/"]
	files = [ "phi0.csv", "phi1.csv"]
	ofname = "aymard-2fol_heatmap.pdf"
	nrow=1
	ncol=2
else:
	print("No relevant case found")
	quit()


# Font sizes and weight
FONT_SIZE_LABEL = 16
FONT_SIZE_TITLE = 18
FONT_SIZE_COLORBAR_LABEL = 16
FONT_SIZE_TICKS = 14
FONT_WEIGHT = "normal"

# ===================




# Load real time values from the summary folder
case_logs = []
all_real_times = []
for folder in dname:
    case_logs.append(pd.read_csv( folder + "case_log.csv"))

for idx in range(len(dname)):
    all_real_times.append(case_logs[idx]["t"].values)
    all_real_times[idx] = np.insert(all_real_times[idx],0, 0.0)


all_mvals = []
for folder in dname:
    all_mvals.append(pd.read_csv(folder + "phi0.csv")["m"].values)

all_matrices = []
epsilon = 1e-3
for folder in dname:
    for fname in files:
        df = pd.read_csv(folder + fname)
        phi_cols = [col for col in df.columns if col.startswith("phi_t")]
        phi_matrix = df[phi_cols].values + epsilon
        all_matrices.append(phi_matrix)


global_min = epsilon
global_max = max(np.max(mat) for mat in all_matrices)
shared_norm = PowerNorm(gamma = 0.3, vmin=global_min, vmax=global_max)

fig, axes = plt.subplots(nrow, ncol, figsize=(10,10), constrained_layout=True)

color_list = ["#3d99cc", "#f2a024", "#74b236", "#eb6235"]
legend_list = ["s0/x0", "s1/x1", "s2/x2", "s3/x3"]

ims = []
for i in range(len(dname)):
    for j in  range(len(files)):
        ax = axes[j,i]
        mvals = all_mvals[i]
        time = all_real_times[i]
        phi_matrix = all_matrices[i*len(files)+j]

        # Create mesh for plotting (centers of rectangles)
        T, M = np.meshgrid(time, mvals)   # shape (N_t, N_m)

        im = ax.pcolormesh(
            T, M, phi_matrix,
            norm=shared_norm,
            shading="auto",
            cmap="binary",
            rasterized=True
        )
        ims.append(im)

        # Logical axes directions
        #ax.set_xlabel("t (days)", fontsize=FONT_SIZE_LABEL, fontweight=FONT_WEIGHT)
        #ax.set_ylabel("m", fontsize=FONT_SIZE_LABEL, fontweight=FONT_WEIGHT)

        if j == 1 and i == 0:   # last row
            ax.set_ylabel("m", fontsize=FONT_SIZE_LABEL, fontweight=FONT_WEIGHT)
        else:
            ax.set_ylabel("")


        if j == len(files) - 1:   # last row
            ax.set_xlabel("t (days)", fontsize=FONT_SIZE_LABEL, fontweight=FONT_WEIGHT)
        else:
            ax.set_xlabel("")


        # Individual subplot title
        ax.set_title(
            rf"$\phi_{{{j}}}$",
            fontsize=FONT_SIZE_TITLE,
            fontweight=FONT_WEIGHT
        )

        ax.tick_params(labelsize=FONT_SIZE_TICKS)
        if i == 1:
            ax.tick_params(axis="y", which="both", left=False, labelleft=False)


        # ====== NEW: add line plot time vs s/x on the same axes ======

        ax2 = ax.twinx()
        t = case_logs[i]["t"].values
        s = case_logs[i][f"s{j}"].values
        x = case_logs[i][f"x{j}"].values
        sx_ratio = s/x

        # Plot s/x vs time on top of heatmap
        ax.plot(
            t,          # x‑axis: time
            sx_ratio,      # y‑axis: s/x
            #color="w",     # white line for visibility
            linewidth=3,
            linestyle="-",
            color=color_list[j],
            label=rf"$s_{j}/x_{j}$"
        )
        ax.legend()
        ax2.tick_params(labelsize=FONT_SIZE_TICKS)
        if j == 1 and i == len(dname) - 1:
            ax2.set_ylabel(r"$s/x$", fontsize=FONT_SIZE_LABEL, fontweight=FONT_WEIGHT)
        else:
            ax2.set_ylabel("")

        if i == 0:
            ax2.tick_params(axis="y", which="both", right=False, labelright=False)




 # One shared colorbar on the right
cbar = fig.colorbar(
    ims[0],
    ax=axes,
    location="right",
    orientation="vertical",
    shrink=0.8
)
cbar.set_label(
    r"$\phi$",
    fontsize=FONT_SIZE_COLORBAR_LABEL,
    weight=FONT_WEIGHT
)

cbar.ax.tick_params(labelsize=FONT_SIZE_TICKS)


# Optional: adjust font globally
plt.rcParams.update({"font.size": FONT_SIZE_LABEL, "font.weight": FONT_WEIGHT})


# save figure
fig.savefig(
    ofname,
    dpi=600,
    bbox_inches="tight",
    format="pdf",
    transparent=True
)
# Png
fig.savefig(ofname.replace(".pdf", ".png"), dpi=600, bbox_inches="tight", transparent=True)

plt.show()

