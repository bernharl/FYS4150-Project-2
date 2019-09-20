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


data = np.loadtxt("eigendata.dat", skiprows=1)
expected = np.array([3, 7, 11, 15])
errors = np.abs(data[:, 1:] - expected)

fig, ax = plt.subplots()
fig.set_size_inches(2 * 2.9, 2 * 1.81134774961)


for i, eigen_error in enumerate(errors.T):
    ax.semilogy(data[:, 0], eigen_error, label=fr"|Err(Eigen={expected[i]})|")

error_lim = np.ones_like(data[:, 0]) * 1e-4
ax.plot(data[:, 0], error_lim, label="Four leading digits")
ax.grid()
ax.set_xlabel(r"$\rho_{max}$")
ax.set_ylabel(r"Absolute error")
fig.tight_layout()
ax.legend()
fig.savefig("../doc/figures/rhoN_plots.eps")
