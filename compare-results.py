from scipy.interpolate import griddata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('prefix')
parser.add_argument('n_samples', nargs=4, type=int)
parser.add_argument('plt_title', type=str)

args = parser.parse_args()

prefix = args.prefix
n_samples = args.n_samples
plt_title = args.plt_title

nipce_df = pd.read_csv(f"data/{prefix}-nipce-{n_samples[0]}-samples/statistics.csv")
mcfew_df = pd.read_csv(f"data/{prefix}-mc-{n_samples[1]}-samples/statistics.csv")
mcmid_df = pd.read_csv(f"data/{prefix}-mc-{n_samples[2]}-samples/statistics.csv")
mcmany_df = pd.read_csv(f"data/{prefix}-mc-{n_samples[3]}-samples/statistics.csv")
det_df = pd.read_csv(f"data/{prefix}-det/result.csv")

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
    plt.plot(mcfew_sample, y, label=f"Monte Carlo ({n_samples[1]} samples)", color="red", linewidth=1)
    plt.plot(mcmid_sample, y, label=f"Monte Carlo ({n_samples[2]} samples)", color='black', linewidth=1)
    plt.scatter(mcmany_sample[::2], y[::2], label=f"Monte Carlo ({n_samples[3]} samples)", c='green', marker='x', s=20)
    plt.plot(nipce_sample, y, label=f"niPCE ({n_samples[0]} samples)", color="blue", linewidth=1)
    
    if varname == "E[u]":
        det_var = det_df["u"].values
        det_sample = griddata(points, det_var, line)
        plt.scatter(det_sample[::4], y[::4], label=f"one deterministic run with mean input", c='m', marker='.', s=20)
    # if varname == "E[u]":
    #     plt.xlim(0.85, 1.05)
    # if varname == "Var[u]":
    #     plt.xlim(0.005, 0.012)
    plt.xlabel(f"{varname} at x={xsample}D")
    plt.ylabel("y/D")
    plt.grid(True)
    plt.legend(loc=2)
    plt.title(plt_title)
    plt.savefig(f"data/{prefix}_sampling_{varname}.png")

plot_var("E[u]")
plot_var("Var[u]")