import numpy as np
import sys
import matplotlib.pyplot as plt
import imageio

def generate_frame(filename):
    unit_length = 1.2e9
    unit_density = 1e-10
    unit_velocity = 1e5
    radii, _, rho, v, P, T, mu, heating, radaccel, heat_eff, n_e, tau_hd, tau_recomb, tau_molec, time, steps = np.loadtxt(filename, unpack=True)
    #plt.figure(0, figsize=(12,8))
    plt.clf()

    fig, axs = plt.subplots(2, 2, sharex='all', num=0)
    axs[0,0].plot(radii, T)
    axs[0,0].set_ylim(100, 7000)
    axs[0,0].set_ylabel("Temperature (K)")

    axs[0,1].plot(radii, v)
    axs[0,1].set_ylim(-10, 30)
    axs[0,1].set_ylabel("Velocity (km/s)")

    axs[1,0].semilogy(radii, rho)
    axs[1,0].set_ylim(1e-10, 10)
    axs[1,0].set_ylabel("Density (1e-10 g/cm^3)")
    axs[1,0].set_xlabel("Radius (Rp)")

    axs[1,1].plot(radii, 4*np.pi*(radii * unit_length)**2* (rho * unit_density) * (v * unit_velocity)/4)
    axs[1,1].set_ylim(-5e10, 5e10)
    axs[1,1].set_ylabel("$\.{M}$ (g/s)")
    axs[1,1].set_xlabel("Radius (Rp)")
    plt.tight_layout()
    plt.savefig("frame.png")

start = int(sys.argv[1])
end = int(sys.argv[2])
if len(sys.argv) == 4:
    skip = int(sys.argv[3])
else:
    skip = 1

writer = imageio.get_writer("movie.mp4", fps=20, mode="I")

for i in range(start, end, skip):
    filename = "data.{}.tab".format(str(i).zfill(4))
    generate_frame(filename)
    
    '''radii, _, rho, v, P, T, mu, heating, radaccel, heat_eff, n_e, tau_hd, tau_recomb, tau_molec, time, steps = np.loadtxt(filename, unpack=True)
    plt.clf()
    plt.plot(radii, T)
    #plt.plot(radii, v)
    plt.ylim(0, 10000)
    #plt.ylim(0, 30)
    plt.xlabel("Radius (Rp)")
    plt.ylabel("Temperature (K)")
    #plt.ylabel("Velocity (km/s)")
    #plt.ylim(0, 2)
    plt.savefig("frame.png")'''
        
    image = imageio.imread("frame.png")
    writer.append_data(image)

writer.close()
