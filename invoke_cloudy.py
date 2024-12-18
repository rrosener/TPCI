import numpy as np
import matplotlib.pyplot as plt
import scipy
import sys

unit_length = 1.2e9
unit_density = 1e-10
Rp = 1.2e9
AMU = 1.67e-24
mu = 0.67

def write(string):
    sys.stdout.write(string + "\n")

filename = "data.{}.tab".format(str(sys.argv[1]).zfill(4))
radii, _, rho, v, P, T, _, heating, radaccel, heat_eff, n_e, tau_hd, tau_recomb, tau_molec, time, steps = np.loadtxt(filename, unpack=True)

plt.loglog(radii, rho)
plt.show()

write('init "exoplanet.ini"')
write('iterate 2')

write('stop thickness {:.4e} linear'.format(10.9 * Rp))
write('dlaw table depth linear')

for i in range(len(radii) - 1, -1, -1):
    depth = radii[-1] - radii[i]
    write("{:.4e} {:.4e}".format(depth * unit_length, rho[i] * unit_density / ( mu * AMU)))
write('end of dlaw')

write('tlaw table depth linear')
for i in range(len(radii) - 1, -1, -1):
    depth = radii[-1] - radii[i]
    write('{:.4e} {:.4e}'.format(depth * unit_length, T[i]))
write('end of tlaw')

'''write('wind advection table depth linear')
for i in range(len(radii)-1, -1, -1):
    depth = radii[-1] - radii[i]
    write('{:.4e} {:.4e}'.format(depth * unit_length, v[i]))
write('end of velocity table')
write('iterate 30')'''
write('iterate 2')
                           
write('element limit off -5')
write('stop temperature linear 5 K')
write('turbulence 1 km/sec no pressure')
write('double optical depth')
write('abundances GASS10')

write('save overview "over.tab"')
write('save pressure "pres.tab" last')
write('save wind "wind.tab"')
write('save dynamics advection "dyna.tab" last')
write('save continuum "continuum.tab" last units Angstrom')
write('save cooling "cool.tab" last')
write('save ages "ages.tab" last')
