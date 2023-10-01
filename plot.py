import numpy as np
import matplotlib.pyplot as plt
import sys

unit_length = 2.15 * 6.378e8
#unit_length = 1.2e9
unit_density = 1e-10
unit_velocity = 1e5

start = int(sys.argv[1])
end = int(sys.argv[2])
if len(sys.argv) == 4:
    skip = int(sys.argv[3])
else:
    skip = 1

for i in range(start, end, skip):
    filename = "data.{}.tab".format(str(i).zfill(4))
    #radii, _, rho, v, P, _,T, mu, heating, radaccel, heat_eff, n_e, tau_hd, tau_recomb, tau_molec, time, steps = np.loadtxt(filename, unpack=True)
    radii, _, rho, v, P, T, mu, heating, radaccel, heat_eff, n_e, tau_hd, tau_recomb, tau_molec, time, steps = np.loadtxt(filename, unpack=True)
    speed_of_sound = np.sqrt(1.38e-16 * T / mu / 1.67e-24)
    n = (1e-10 * rho) / (2.3 * 1.67e-24)
    sigma = 2e-15
    mfp = 1.0 / (n * sigma)
    #plt.semilogy(radii, mfp)
    #plt.figure()
    #plt.plot(radii, T)
    plt.figure(0)
    plt.plot(radii, v)
    plt.title("Velocity (km/s)")
    plt.figure(1)
    plt.semilogy(radii, rho)    
    plt.title("Density (g/cm^3)")

    plt.figure(2)
    plt.semilogy(radii, P)
    plt.title("Pressure (dyn/cm^2)")

    plt.figure(3)
    plt.plot(radii, T)
    plt.title("Temperature (K)")
    #plt.title(str(i))

    plt.figure(4)
    plt.plot(radii, 4*np.pi*(radii * unit_length)**2* (rho * unit_density) * (v * unit_velocity))
    plt.ylim(0, 1e10)
    plt.title("Mass loss rate (g/s)")
    
   # plt.figure(5)
   # plt.plot(radii, P**(1-5./3) * T**(5./3))

    #dln_rho_dr = np.diff(np.log(rho))
    #dln_rho_dr = 


plt.show()
