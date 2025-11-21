from scipy.interpolate import griddata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df10 = pd.read_csv("data/simple_res10.csv")
df20 = pd.read_csv("data/simple_res20.csv")
df50 = pd.read_csv("data/simple_res50.csv")

pt10 = df10[['x', 'y']].values
pt20 = df20[['x', 'y']].values
pt50 = df50[['x', 'y']].values

xsample = 10
y = np.linspace(df10['y'].min(), df10['y'].max(), 100)
x = np.full_like(y, xsample)
line = np.column_stack((x, y))

u10 = df10['u'].values
u20 = df20['u'].values
u50 = df50['u'].values

prof10 = griddata(pt10, u10, line)
prof20 = griddata(pt20, u20, line)
prof50 = griddata(pt50, u50, line)

plt.figure(figsize=(5,6))
plt.scatter(prof50[::5], y[::5], label='dx=D/50', c='red', marker='+', s=50)
plt.scatter(prof20[::3], y[::3], label='dx=D/20', c='green', marker='.', s=60)
plt.plot(prof10, y, label='dx=D/10', color='blue', linewidth=1)
plt.xlabel("u")
plt.ylabel("y/D")
plt.title("velocity profile at x=10D")
plt.legend()
plt.grid(True)
plt.savefig("data/res-sensitivity.png")