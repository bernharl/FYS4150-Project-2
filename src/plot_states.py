import numpy as np 
import matplotlib.pyplot as plt


fonts = {
        "font.family": "serif",
        "axes.labelsize": 10,
        "font.size": 10,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
}
plt.rcParams.update(fonts)
omega_r = np.loadtxt("omegas.txt", skiprows = 1)
rhoN = np.loadtxt("omegas.txt", max_rows = 1)
eigenvecs = np.loadtxt("ground_states.txt", skiprows=2)
eigvals = np.loadtxt("lambdas.txt", skiprows=2)

n = len(eigenvecs[:,0]) 
h = rhoN / (1 + len(eigenvecs[:,0]))
rho = np.linspace(h, n*h, len(eigenvecs[:,0]))

fig, ax = plt.subplots()
ax.plot(rho, np.abs(eigenvecs))
plt.show()



