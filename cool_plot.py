import numpy as np
import matplotlib.pyplot as plt
import sys

max_R = 14.88
R_earth = 6.378136e8
Rp = 2.9 * R_earth

radii = []
depths = []
cooling = {}
heating_rate = []
cooling_rate = []

for line in open(sys.argv[1]):
    if line[0] == '#': continue
    elements = line.split("\t")
    depths.append(float(elements[0]))
    radii.append(14.88 * Rp - depths[-1])
    heating_rate.append(float(elements[2]))
    cooling_rate.append(float(elements[3]))

    pos = 4
    while pos < len(elements):
        name = elements[pos]
        if name not in cooling:
            cooling[name] = []
        cooling[name].append([radii[-1], float(elements[pos + 1])])
        pos += 2

radii = np.array(radii)        
#radii = 14.88 * Rp - np.array(depths)
#print(cooling)
colors = {"ISrcolH H  0.0": '#1f77b4', "FF c 0.0": '#ff7f0e', "ISrcolHeHe 1.0": '#2ca02c', "Fe 2 0.0": '#d62728', "Mg 2 2802.7": '#9467bd', "ISclinH H  0.0": '#8c564b', "Fe 6 0.0": '#e377c2', "O  3 5007.0": '#7f7f7f', "Si 3 1892.0": '#bcbd22', "Mg 2 2795.5": '#17becf', "Ca 2 7306.0": "k", "Na 1 5891.9": "b", "Fe 1 240359.5": "pink", "ISrcolH He 1.0": "gray"}
replacements = {"ISrcolH H  0.0": "H II recomb",
                "FF c 0.0": "free-free",
                "ISrcolHeHe 1.0": "He II recomb",
                "ISclinH H  0.0": "H line",
                "ISrcolH He 1.0": "He III recomb",
                "Fe 2 0.0": "Fe 2",
                "Fe 6 0.0": "Fe 6"
}

for name in cooling:
    cooling[name] = np.array(cooling[name])
    if np.max(cooling[name][:,1]) > 0.15:
        label = name
        if name in replacements:
            label = replacements[name]
        plt.semilogx(cooling[name][:,0] / Rp, cooling[name][:,1], label=label, color=colors[name])
        
plt.legend(frameon=False)
plt.xlabel("Radius (Rp)", fontsize=12)
plt.ylabel("Cooling contribution", fontsize=12)
plt.xlim(1, 15)
plt.ylim(0, 0.8)

'''plt.figure()
plt.loglog(radii / Rp, heating_rate, label="Rad heating")
plt.loglog(radii / Rp, cooling_rate, label="Rad cooling")

cooling_radii = np.load("adiabatic_cooling_radii.npy")[1:-1]
ad_cooling = np.load("adiabatic_cooling_0.npy")[1:-1]
advec_cooling = np.load("advec_cooling_0.npy")[1:-1]
plt.loglog(cooling_radii / Rp, ad_cooling, label="Adiabatic cooling")
plt.loglog(cooling_radii / Rp, np.abs(advec_cooling), label="|Advection|")
plt.xlim(1, 15)
plt.ylim(1e-11, 2e-4)
plt.xlabel("Radius (Rp)", fontsize=12)
plt.ylabel("Rate (erg s$^{-1}$ cm$^{-3}$)", fontsize=12)
plt.title("Z=1")
plt.legend()
plt.tight_layout()
plt.savefig("thermal_terms_1Z.png")'''
plt.show()    
