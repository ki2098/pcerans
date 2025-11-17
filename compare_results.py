from scipy.interpolate import griddata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

prefix = sys.argv[1]


nipce_df = pd.read_csv(f"data/{prefix}_nipce_sample.statistics.csv")
mcfew_df = pd.read_csv(f"data/{prefix}_mc_sample.statistics-few.csv")
mcmid_df = pd.read_csv(f"data/{prefix}_mc_sample.statistics-mid.csv")
mcmany_df = pd.read_csv(f"data/{prefix}_mc_sample.statistics-many.csv")
det_df = pd.read_csv(f"data/{prefix}_det.csv")

points = nipce_df[["x", "y"]].values

xsample = 10
y = np.linspace(nipce_df["y"].min(), nipce_df["y"].max(), 100)
x = np.full_like(y, xsample)
line = np.column_stack((x, y))

def plot_var(varname:str):
    nipce_var = nipce_df[varname].values
    nipce_sample = griddata(points, nipce_var, line)

    mcfew_var = mcfew_df[varname].values
    mcfew_sample = griddata(points, mcfew_var, line)

    mcmid_var = mcmid_df[varname].values
    mcmid_sample = griddata(points, mcmid_var, line)

    mcmany_var = mcmany_df[varname].values
    mcmany_sample = griddata(points, mcmany_var, line)

    plt.figure(figsize=(5,6))
    plt.plot(mcfew_sample, y, label=f"Monte Carlo ({sys.argv[2]} samples)", color="red", linewidth=1)
    plt.plot(mcmid_sample, y, label=f"Monte Carlo ({sys.argv[3]} samples)", color='black', linewidth=1)
    plt.scatter(mcmany_sample[::2], y[::2], label=f"Monte Carlo ({sys.argv[4]} samples)", c='green', marker='x', s=20)
    plt.plot(nipce_sample, y, label=f"niPCE ({sys.argv[5]} samples)", color="blue", linewidth=1)
    
    # if varname == "E[u]":
    #     det_var = det_df["u"].values
    #     det_sample = griddata(points, det_var, line)
    #     plt.scatter(det_sample[::4], y[::4], label=f"deterministic run uin=E[uin]", c='black', marker='x', s=20)
    # if varname == "E[u]":
    #     plt.xlim(0.85, 1.05)
    # if varname == "Var[u]":
    #     plt.xlim(0.005, 0.012)
    plt.xlabel(f"{varname} at x={xsample}")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend(loc=2)
    plt.title(sys.argv[6])
    plt.savefig(f"data/{prefix}_sampling_{varname}.png")

plot_var("E[u]")
plot_var("Var[u]")