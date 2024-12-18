import numpy as np
import matplotlib
import matplotlib.pyplot as plt

e = 4.8032e-10
m_e = 9.11e-28
c = 3e10
A = 1.0216e7
m_He = 4 * 1.67e-24
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


def plot(radii, Te, wind_vel, hden, n1, n3, HI, HeI, label, axs): 
    axs[0,0].plot(radii/Rp, Te, label=label)
    axs[0,0].set_ylabel("Temperature (K)")
    axs[0,0].set_ylim(0, 15000)
    
    axs[0,1].plot(radii/Rp, wind_vel/1e5)
    axs[0,1].set_ylabel("Velocity (km/s)")
    axs[0,1].set_ylim(0, 30)
    
    axs[1,0].plot(radii/Rp, hden)
    axs[1,0].set_ylabel("$n_H$ (cm$^{-3}$)")
    axs[1,0].set_yscale("log")
    
    axs[1,1].plot(radii/Rp, n3)
    axs[1,1].set_ylabel("$n_{He*}$ (cm$^{-3}$)")
    axs[1,1].set_ylim(1, 1e5)
    axs[1,1].set_yscale("log")

    axs[2,0].plot(radii/Rp, HI)
    axs[2,0].set_ylabel("Neutral H frac.")
    axs[2,0].set_xlabel("Radius (R$_p$)")
    
    axs[2,1].plot(radii/Rp, HeI)
    axs[2,1].set_ylabel("Neutral He frac.")
    axs[2,1].set_xlabel("Radius (R$_p$)")

    #axs[2,1].plot(radii/Rp, n3 / (n1 + n3))
    #axs[2,1].set_ylabel("He3/He")
    #axs[2,1].set_xlabel("Radius (R$_p$)")
    #axs[2,1].set_yscale("log")

    plt.xlim(0.99,10)

    for ax in axs.flat:
        ax.set_xscale("log")
        ax.set_xticks([])
        ax.set_xticks(ticks=[1, 2, 3, 4, 5, 10])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.FixedFormatter(["1", "2", "3", "4", "5", "10"]))
        ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

def plot_quantities(dirname, number, axs, label, multiple):
    radii = []
    n1 = []
    n3 = []
    ne = []

    for line in open("{}/cl_data.{}.pops.tab".format(dirname, number)):
        elements = line.split()
        if len(elements) > 10 and elements[1] == "He":
            radii.append(float(elements[0]))
            n1.append(float(elements[2]))
            n3.append(float(elements[3]))
        if len(elements) == 3 and elements[1] == "e-":
           ne.append(float(elements[2])) 

    if np.max(radii) / Rp + 1 <= 12:
        max_R = 11.92451
    elif np.max(radii) / Rp + 1 > 14:
        max_R = 14.9 #88
    else:
        assert(False)

    radii = (max_R*Rp - np.array(radii))[::-1]
    n1 = np.array(n1)[::-1]
    n3 = np.array(n3)[::-1]
    
    wind_radii, wind_vel = np.loadtxt("{}/cl_data.{}.wind.tab".format(dirname, number), unpack=True, usecols=(1,2))
    wind_radii = max_R*Rp - wind_radii[::-1]
    wind_radii = wind_radii[:len(radii)]
    wind_vel = -1 * wind_vel[::-1][:len(radii)]

    data = np.loadtxt("{}/cl_data.{}.over.tab".format(dirname, number))[::-1][0:len(radii)]
    _, Te, Htot, hden, eden, H_molec_ratio, HI, HII, HeI, HeII, HeIII = data.T[0:11]
    plot(radii, 10**Te, wind_vel, 10**hden, n1, n3, 10**HI, 10**HeI, label, axs)

        
fig, axs = plt.subplots(3, 2, sharex='all', figsize=(9,7))
plot_quantities("55cnc_final", "0123", axs, "90% H, 10% He", 1./9)
plot_quantities("He_final", "0123", axs, "1% H, 99% He", 99)

#Handle Salz et al
radii, rho, vel, P, T, mu, heating, heating_frac, HI, HII, HeI, HeII, HeIII = np.loadtxt("55cnc/salz_55cnc.txt", unpack=True)

plot(radii*Rp, T, vel*1e5, rho/1.67e-24*0.9, np.zeros(len(radii)), np.zeros(len(radii)), HI, HeI, "Salz et al 2016", axs)

#Handle 2.5D
radii, rho, T, v, n_H, n_H_ion, n_H2, n_He, n_He_ion, n3 = np.loadtxt("profile_fiducial.txt", unpack=True)
plot(radii * R_earth, T, v, n_H, n_He, n3, n_H / (n_H_ion + n_H), n_He / (n_He + n_He_ion), "2.5D", axs)

handles, labels = axs[0,0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=4)
#plt.tight_layout()
plt.show()
