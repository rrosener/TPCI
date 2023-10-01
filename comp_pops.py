import numpy as np
import matplotlib
import matplotlib.pyplot as plt

e = 4.8032e-10

m_e = 9.11e-28
c = 3e10
A = 1.0216e7
m_H = 1.67e-24
m_He = 4 * m_H
k_B = 1.38e-16
R_sun = 7e10
M_jup = 1.898e30
M_earth = 5.97e27
R_jup = 7.1e9
R_earth = 6.378e8
parsec = 3.086e18
AU = 1.496e13
Rp = 1.3e9
Rs = 0.943 * R_sun


def plot(radii, Te, wind_vel, n_H, n_He, n_e, label, axs):
    rhos = n_H * m_H + n_He * m_He
    mus = rhos / (n_H + n_He + n_e) / m_H
    
    axs[0,0].plot(radii/Rp, Te / mus, label=label)
    axs[0,0].set_ylabel("$T/\mu$ (K/amu)")
    #axs[0,0].set_ylim(0, 15000)
    
    axs[0,1].plot(radii/Rp, mus)
    axs[0,1].set_ylabel("$\mu$ (amu)")
    
    axs[1,0].plot(radii/Rp, n_H)
    axs[1,0].set_ylabel("$n_H$ (cm$^{-3}$)")
    axs[1,0].set_yscale("log")
    
    axs[1,1].plot(radii/Rp, n_He)
    axs[1,1].set_ylabel("$n_{He}$ (cm$^{-3}$)")
    #axs[1,1].set_ylim(1, 1e5)
    axs[1,1].set_yscale("log")

    axs[2,0].plot(radii/Rp, n_e)
    axs[2,0].set_ylabel("$n_e$ (cm$^{-3}$")
    axs[2,0].set_xlabel("Radius (R$_p$)")
    axs[2,0].set_yscale("log")
    
    axs[2,1].plot(radii/Rp, rhos)
    axs[2,1].set_ylabel("Density (g/cm^3)")
    axs[2,1].set_xlabel("Radius (R$_p$)")
    axs[2,1].set_yscale("log")

    plt.xlim(0.99,10)

    for ax in axs.flat:
        ax.set_xscale("log")
        ax.set_xticks([])
        ax.set_xticks(ticks=[1, 2, 3, 4, 5, 10])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.FixedFormatter(["1", "2", "3", "4", "5", "10"]))
        ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

def plot_quantities(dirname, number, axs, label):
    radii = []
    n_HeI = []
    n_HeII = []
    n_HeIII = []
    n_e = []
    n3 = []

    n_HI = []
    n_HII = []
    
    for line in open("{}/cl_data.{}.pops.tab".format(dirname, number)):
        elements = line.split()
        if len(elements) > 10 and elements[1] == "He":
            radii.append(float(elements[0]))
            n_HeI.append(np.sum([float(f) for f in elements[2:]]))
            n3.append(float(elements[3]))
        if len(elements) > 2 and elements[1] == "He+":
            n_HeII.append(np.sum([float(f) for f in elements[2:]]))
        if len(elements) > 2 and elements[1] == "He+2":
            n_HeIII.append(np.sum([float(f) for f in elements[2:]]))
        if len(elements) == 3 and elements[1] == "e-":
           n_e.append(float(elements[2]))
        if len(elements) > 2 and elements[1] == "H":
            n_HI.append(np.sum([float(f) for f in elements[2:]]))
        if len(elements) > 2 and elements[1] == "H+":
            n_HII.append(np.sum([float(f) for f in elements[2:]]))
    
    if np.max(radii) / Rp + 1 <= 12:
        max_R = 11.92451
    elif np.max(radii) / Rp + 1 > 14:
        max_R = 14.9 #88
    else:
        assert(False)

    radii = (max_R*Rp - np.array(radii))#[::-1]
    
    wind_radii, wind_vel = np.loadtxt("{}/cl_data.{}.wind.tab".format(dirname, number), unpack=True, usecols=(1,2))
    #wind_radii = max_R*Rp - wind_radii[::-1]
    #wind_radii = wind_radii[:len(radii)]
    wind_vel = -1 * wind_vel[-len(radii):]

    data = np.loadtxt("{}/cl_data.{}.over.tab".format(dirname, number))[-len(radii):]
    _, Te, Htot, hden, eden, H_molec_ratio, HI, HII, HeI, HeII, HeIII = data.T[0:11]

    plot(radii, 10**Te, wind_vel,
         np.array(n_HI) + np.array(n_HII),
         np.array(n_HeI) + np.array(n_HeII) + np.array(n_HeIII),
         np.array(n_e),
         label, axs)

        
fig, axs = plt.subplots(3, 2, sharex='all', figsize=(9,7))
plot_quantities("55cnc_final", "0123", axs, "90% H, 10% He")
plot_quantities("He_final", "0123", axs, "1% H, 99% He")


handles, labels = axs[0,0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=4)
#plt.tight_layout()
plt.show()
