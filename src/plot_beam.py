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

analytical_data = np.loadtxt("bbeam_analytical.dat", skiprows=1) 
numerical_data = np.loadtxt("bbeam_num.dat", skiprows=1)
N = numerical_data[:,0] 

rel_err1 = np.abs((numerical_data[:,1] - analytical_data[:,1])
                   / analytical_data[:,1])
rel_err2 = np.abs((numerical_data[:,2] - analytical_data[:,2])
                   / analytical_data[:,2])
rel_err3 = np.abs((numerical_data[:,3] - analytical_data[:,3])
                   / analytical_data[:,3])

labels = [rf"$\lambda_{i+1}$" for i in range(3)]
fig, ax = plt.subplots()
fig.set_size_inches(2 * 2.9, 2 * 1.81134774961)
ax.semilogy(N, rel_err1)
ax.semilogy(N, rel_err2)
ax.semilogy(N, rel_err3)
ax.set_xlabel(r"$N$")
ax.set_ylabel(r"Relative error log($\vert \frac{\lambda_{num}-\lambda_{analy}}{\lambda_{analy}}\vert$)")
ax.grid()
ax.legend(labels)
fig.tight_layout()
fig.savefig("../doc/figures/eigenrelerr.eps", dpi=1000)
