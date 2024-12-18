import numpy as np
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

def plot_pops(dirname, number, color, label):
    radii = []
    n_HeI = []
    n_HeII = []
    n_HeIII = []
    ne = []
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
           ne.append(float(elements[2]))
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

    radii = (max_R*Rp - np.array(radii))


    
    
    n_HeI = np.array(n_HeI)
    n_HeII = np.array(n_HeII)
    n_HeIII = np.array(n_HeIII)
    n_ion = n_HeII + n_HeIII
    ne = np.array(ne)

    alpha3 = 3.4e-13
    q31a = 2.6e-8
    ion_rate3 = 3.5
    
    n3_pred_deexcite = n_ion * alpha3 / q31a
    n3_pred_ion = n_ion * ne * alpha3 / ion_rate3
    n3_pred = (n_ion * ne * alpha3) / (ion_rate3 + ne * q31a)
    print(len(radii), len(n_HeI), len(n_HeII), len(n_HeIII), len(ne))

    #plt.loglog(radii/Rp, n_HeII + n_HeIII)
    #plt.loglog(radii / Rp, (n_HeII + n_HeIII) * ne * 1e-10)
    plt.loglog(radii / Rp, n3, color=color, label=label + " (TPCI)")
    plt.loglog(radii / Rp, n3_pred, color=color, linestyle='--', label=label + " (analytic)")
    plt.xlabel("Radius (Rp)")
    plt.ylabel("n_He* (cm^-3)")
    #plt.loglog(radii / Rp, n3, color=color, linestyle="--")
    #plt.show()


#dirname = "55cnc_final"
#number = "0125"    
plot_pops("55cnc_final", "0123", "b", "H-dom")
plot_pops("He_final", "0123", "orange", "He-dom")
plt.legend()
plt.show()
