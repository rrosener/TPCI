import scipy.special
import numpy as np
import scipy.integrate
import scipy.interpolate
import matplotlib.pyplot as plt
import sys
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import interp1d

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
R_earth = 6.378136e8
parsec = 3.086e18
AU = 1.496e13
#Rp = 1.3e9
#Rs = 0.943 * R_sun
#T0 = 5000

#Rp = 4.36 * R_earth
#Rs = 0.683 * R_sun

#Rp = 4.57 * R_earth
#Rs = 0.547 * R_sun
#T0 = 7000

#HD 63433 system
#Rp = 2.67 * R_earth
#Rp = 2.15 * R_earth
#Rs = 0.912 * R_sun
#T0 = 8000

#TOI 560
#Rp = 2.9 * R_earth
#Rs = 0.665 * R_sun
#T0 = 7000
#hill_radius = 11 #planetary radii

#HD 189733
Rp = 12.7 * R_earth
Rs = 0.78 * R_sun
T0 = 7000
hill_radius = 6

def get_tau(wavenum, b, n3_interp, T0, max_r, v_interp, epsilon=0.01):
    line_wavenums, fs = np.loadtxt("../He_line_data", unpack=True)
    assert(len(line_wavenums) == 3)
    #wavenums = 9230.868568
    #fs = 1.7974e-1
    
    sigma_0 = np.pi * e**2 * fs / m_e / c**2
    doppler_broadening = np.sqrt(k_B * T0 / m_He) * line_wavenums / c
    
    def integrand(theta, line_index):
        r = b / np.cos(theta)
        dtheta_to_dr = np.abs(b*np.sin(theta)/np.cos(theta)**2)
        gamma = A / 4 / np.pi / c
        lorentz = 1.0/np.pi * gamma / ((wavenum - line_wavenums[line_index])**2 + gamma**2)
        v_offset = v_interp(r) * np.sin(theta) 
        
        z = (wavenum - line_wavenums[line_index] - line_wavenums[line_index]/c*v_offset + gamma * 1j) / doppler_broadening[line_index] / np.sqrt(2)
        profile = np.real(scipy.special.wofz(z)) / doppler_broadening[line_index] / np.sqrt(2*np.pi)
        per_line = n3_interp(r) * sigma_0[line_index] * profile * r / np.sqrt(r**2 - b**2) * dtheta_to_dr
        return per_line

    max_theta = np.arccos(b / max_r)
    min_theta = -max_theta

    #Rough integration
    thetas = np.linspace(min_theta + 1e-2, max_theta - 1e-2, 250)
    integrands = np.sum([integrand(thetas, i) for i in range(len(line_wavenums))], axis=0)
    result = np.trapz(integrands, thetas)
    return result


def predict_depths(wavenums, r_mesh, v_interp, n3_interp, T0, max_r_over_Rp=9.9):
    radii = np.linspace(Rp, max_r_over_Rp * Rp, 200)
    dr = np.median(np.diff(radii))

    tot_extra_depth = 0
    transit_spectrum = []

    for wavenum in wavenums:
        tot_extra_depth = 0
        for r in radii:
            tau = get_tau(wavenum, r, n3_interp, T0, np.max(radii), v_interp)
            extra_depth = 2 / Rs**2 * r * dr * (1 - np.exp(-tau))
            tot_extra_depth += extra_depth        

        #print(wavenum, tot_extra_depth * 1e6)
        transit_spectrum.append(tot_extra_depth)

    return np.array(transit_spectrum)

#number = sys.argv[1]
labels = ["Z=1", "Z=10", "Z=30", "Z=100"]
colors = ["k", "#5e110f", "#ad201b", "#ff2f28"]
for i, number in enumerate(sys.argv[1:]):
    radii = []
    n1 = []
    n3 = []
    ne = []

    for line in open("cl_data.{}.pops.tab".format(number)):
        elements = line.split()
        if len(elements) > 10 and elements[1] == "He":
            radii.append(float(elements[0]))
            n1.append(float(elements[2]))
            n3.append(float(elements[3]))
        if len(elements) == 3 and elements[1] == "e-":
            ne.append(float(elements[2]))

    #print(radii)
    if np.max(radii) / Rp + 1 <= 12:
        #max_R = 11.92451
        max_R = np.max(radii) / Rp + 1
        #print(max_R)
    elif np.max(radii) / Rp + 1 > 14:
        max_R = 14.88
    else:
        assert(False)

    adjust_factor = 2 #for 3D structure    
    radii = max_R*Rp - np.array(radii)
    n1 = np.array(n1) / adjust_factor #Account for 3D structure 
    n3 = np.array(n3) / adjust_factor
    ne = np.array(ne) / adjust_factor
    #n3 /= 15**2

    radii = radii[::-1]
    n1 = n1[::-1]
    n3 = n3[::-1]
    ne = ne[::-1]

    #np.save("tpci_radii.npy", radii)
    #np.save("tpci_n3.npy", n3)

    wind_radii, wind_vel = np.loadtxt("cl_data.{}.wind.tab".format(number), unpack=True, usecols=(1,2))
    wind_radii = max_R*Rp - wind_radii
    wind_radii = wind_radii[::-1][:len(radii)]
    wind_vel = -1 * wind_vel[::-1][:len(radii)] / 2 #Account for 3D structure
    #plt.plot(wind_radii/Rp, wind_vel/1e5)
    #plt.show()

    assert(np.allclose(radii, wind_radii))
    #print(wind_radii)
    #print(radii)

    res = 375000
    wavenums = np.exp(np.arange(np.log(9230), np.log(9233), 1./res))
    depths = predict_depths(
        wavenums, radii,
        interp1d(radii, wind_vel, bounds_error=False, fill_value=(0, wind_vel[-1])),
        interp1d(radii, n3, bounds_error=False, fill_value=(n3[0], n3[-1])),
        T0, max_r_over_Rp=hill_radius)
    print(labels[i], 1e2 * np.max(depths), -1e3 * np.trapz(depths, 1e8/wavenums))

    plt.figure(0)
    #plt.loglog(radii/Rp, n1)
    plt.loglog(radii/Rp, n3, label=labels[i], color=colors[i])
    plt.xlabel("Radius (Rp)", fontsize=12)
    plt.ylabel("$n_{He^*}$ (cm$^{-3}$)", fontsize=12)
    plt.xlim(1, 15)
    plt.ylim(1e-1, 1e4)
    plt.tight_layout()
    
    plt.figure(1)
    #plt.plot(1e8/wavenums, 1e2 * depths, label=labels[i], color=colors[i])
    plt.plot(1e8/wavenums, 1e2 * gaussian_filter(depths, 375000/37500/2.355), label=labels[i], color='r')
    plt.xlabel("Wavelength ($\AA$)", fontsize=12)
    plt.ylabel("Excess absorption (%)", fontsize=12)
    plt.tight_layout()

    plt.figure(2)
    plt.loglog(radii/Rp, ne, label=labels[i], color=colors[i])
    plt.xlabel("Radius (Rp)", fontsize=12)
    plt.ylabel("$n_e$ (cm$^{-3}$)", fontsize=12)
    plt.xlim(1, 15)
    plt.tight_layout()
    #np.save("depths.npy", depths)

plt.figure(0)
plt.legend()
plt.figure(1)
plt.legend()
plt.show()

#7.5 Mearth: 0105, 0144, 0195, 0210
#9 Mearth: 0170, 0213, 1009, 1170
