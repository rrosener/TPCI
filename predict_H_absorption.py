import scipy.special
import numpy as np
import scipy.integrate
import scipy.interpolate
import matplotlib.pyplot as plt
import sys
import scipy.ndimage.filters
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
R_earth = 6.378e8
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

Rp = 2.67 * R_earth
#Rp = 2.15 * R_earth
Rs = 0.912 * R_sun
T0 = 10000

#TOI 560
#Rp = 2.8 * R_earth
#Rs = 0.665 * R_sun
#T0 = 7000
#hill_radius = 11

def get_tau(wavenum, b, n3_interp, T0, max_r, v_interp, epsilon=0.01):
    line_wavenums, fs = np.loadtxt("../H_line_data", unpack=True)
    line_wavenums = np.atleast_1d(line_wavenums)
    fs = np.atleast_1d(line_wavenums)
    #assert(len(line_wavenums) == 1)
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
        assert(np.all(per_line >= 0))
        return per_line

    max_theta = np.arccos(b / max_r)
    min_theta = -max_theta
    assert(max_theta - 1e-2 > min_theta + 1e-2)

    #Rough integration
    thetas = np.linspace(min_theta + 1e-2, max_theta - 1e-2, 250)
    integrands = np.sum([integrand(thetas, i) for i in range(len(line_wavenums))], axis=0)
    result = np.trapz(integrands, thetas)
    if result < 0:
        assert(False)
    assert(result >= 0)
    return result


def predict_depths(wavenums, r_mesh, v_interp, n3_interp, T0, max_r_over_Rp=9.9):
    radii = np.linspace(Rp, max_r_over_Rp * Rp, 200)
    dr = np.median(np.diff(radii))

    tot_extra_depth = 0
    transit_spectrum = []

    for wavenum in wavenums:
        tot_extra_depth = 0
        for r in radii[1:-1]:
            tau = get_tau(wavenum, r, n3_interp, T0, np.max(radii), v_interp)
            extra_depth = 2 / Rs**2 * r * dr * (1 - np.exp(-tau))
            assert(extra_depth >= 0)
            tot_extra_depth += extra_depth        

        #print(wavenum, tot_extra_depth * 1e6)
        transit_spectrum.append(tot_extra_depth)

    return np.array(transit_spectrum)

number = sys.argv[1]
radii = []
n_H = []
n_HII = []

for line in open("cl_data.{}.pops.tab".format(number)):
    elements = line.split()
    if len(elements) > 10 and elements[1] == "H":
        radii.append(float(elements[0]))
        n_H.append(float(elements[2]))
    if len(elements) == 3 and elements[1] == "H+":
        n_HII.append(float(elements[2]))

if np.max(radii) / Rp + 1 <= 12:
    max_R = 11.92451
elif np.max(radii) / Rp + 1 > 14:
    max_R = 14.88
else:
    assert(False)
    
radii = max_R*Rp - np.array(radii)
n_H = np.array(n_H) / 2 #Account for 3D structure
n_HII = np.array(n_HII) / 2
#n_H /= 15**2

radii = radii[::-1]
n_H = n_H[::-1]
n_HII = np.array(n_HII)[::-1]

wind_radii, wind_vel = np.loadtxt("cl_data.{}.wind.tab".format(number), unpack=True, usecols=(1,2))
wind_radii = max_R*Rp - wind_radii
wind_radii = wind_radii[::-1][:len(radii)]
wind_vel = -1 * wind_vel[::-1][:len(radii)] / 2 #Account for 3D structure
#plt.plot(wind_radii/Rp, wind_vel/1e5)
#plt.show()

assert(np.allclose(radii, wind_radii))
#print(wind_radii)
#print(radii)

recomb_rates = n_HII**2 * 2.58e-13 #at 10,000 K
advec_rates = -wind_vel[1:] * np.diff(n_H) / np.diff(radii)
advec_rates = np.append([0], advec_rates)
plt.semilogy(radii / Rp, recomb_rates, label="Recomb")
plt.semilogy(radii / Rp, advec_rates, label="Advec")
plt.legend()
plt.xlabel("Radius (planetary radii)")
plt.ylabel("Rates (cm^-6 s^-1)")
plt.show()


res = 37500
wavenums = np.exp(np.arange(np.log(80000), np.log(84000), 1./res))
depths = predict_depths(
    wavenums, radii,
    interp1d(radii, wind_vel, bounds_error=False, fill_value=(0, wind_vel[-1])),
    interp1d(radii, n_H, bounds_error=False, fill_value=(n_H[0], n_H[-1])),
    T0, max_r_over_Rp=hill_radius)

plt.semilogy(radii/Rp, n_H)
plt.figure()
plt.plot(1e4/wavenums, depths)
plt.ylim(0, 1)
plt.show()
