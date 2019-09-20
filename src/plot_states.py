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


n = len(eigenvecs[:,0])
h = rhoN / (1 + len(eigenvecs[:,0]))
rho = np.linspace(h, n*h, len(eigenvecs[:,0]))
labels = [rf"$\omega=${i}" for i in omega_r]
fig, ax = plt.subplots()
fig.set_size_inches(2 * 2.9, 2 * 1.81134774961)
ax.plot(rho, np.abs(eigenvecs))
ax.set_xlabel(r"$\rho$")
ax.set_ylabel(r"Ground state eigenvector $\psi$")
ax.grid()
ax.legend(labels)
fig.tight_layout()
fig.savefig("../doc/figures/eigenstates.eps", dpi=1000)
