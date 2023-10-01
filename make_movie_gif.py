import numpy as np
import sys
import matplotlib.pyplot as plt
import imageio

def generate_frame(filename):
    unit_length = 1.2e9
    unit_density = 1e-10
    unit_velocity = 1e5
    radii, _, rho, v, P, T, mu, heating, radaccel, heat_eff, n_e, tau_hd, tau_recomb, tau_molec, time, steps = np.loadtxt(filename, unpack=True)
    plt.figure(0, figsize=(5,4))
    plt.clf()    
    plt.plot(radii, v)
    plt.ylim(-5, 40)
    plt.xlabel("Radius (Rp)")
    plt.ylabel("Velocity (km/s)")

    plt.savefig("frame.png")

start = int(sys.argv[1])
end = int(sys.argv[2])
if len(sys.argv) == 4:
    skip = int(sys.argv[3])
else:
    skip = 1

writer = imageio.get_writer("poster_movie.gif", fps=20, mode="I")

for i in range(start, end, skip):
    filename = "data.{}.tab".format(str(i).zfill(4))
    generate_frame(filename)
    
    image = imageio.imread("frame.png")
    writer.append_data(image)

writer.close()
