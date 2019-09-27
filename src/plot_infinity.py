import numpy as np
import matplotlib.pyplot as plt

# Setting fonts for pretty plot
fonts = {
    "font.family": "serif",
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
}
plt.rcParams.update(fonts)

# Reading data, defining errors
data = np.loadtxt("eigendata.dat", skiprows=1)
expected = np.array([3, 7, 11, 15])
errors = np.abs(data[:, 1:] - expected)

fig, ax = plt.subplots()
# The golden ratio for nice plots size, also based on width of tex document
fig.set_size_inches(2 * 2.9, 2 * 1.81134774961)

# Plotting dats
for i, eigen_error in enumerate(errors.T):
    ax.semilogy(data[:, 0], eigen_error, label=fr"|Err($\lambda$={expected[i]})|")

# We want errors below 1e-3
error_lim = np.ones_like(data[:, 0]) * 1e-3
ax.plot(data[:, 0], error_lim, label="Four leading digits")
# Finding best RhoN for lambda=15
min_pos = np.argmin(errors[:, -1])
# ax.scatter(data[min_pos, 0], errors[min_pos, -1], label="Best")
ax.plot(
    (data[min_pos, 0], data[min_pos, 0]),
    (0, np.max(errors[:, -1])),
    "--",
    label=rf"Best $\rho_{{max}} \approx {data[min_pos, 0]:.2f}$",
)
ax.grid()
ax.set_xlabel(r"$\rho_{max}$")
ax.set_ylabel(r"Absolute error")
fig.tight_layout()
ax.legend(loc=1)
# Saving high quality figure
fig.savefig("../doc/figures/rhoN_plots.eps", dpi=1000)
