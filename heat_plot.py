import numpy as np
import matplotlib.pyplot as plt
import sys

max_R = 14.88
R_earth = 6.378136e8
#Rp = 2.9 * R_earth
Rp = 4.36 * R_earth

radii = []
depths = []
cooling = {}
heating_rate = []
cooling_rate = []

for line in open(sys.argv[1]):
    if line[0] == '#': continue
    if "LHeat" in line: continue
    elements = line.split("\t")
    print(elements)
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

heating_rate = np.array(heating_rate)
cooling_rate = np.array(cooling_rate)

#radii = 14.88 * Rp - np.array(depths)
#print(cooling)
colors = {"He 2": '#1f77b4', "H  1": '#ff7f0e', "He 1": '#2ca02c', "Hlin": '#d62728', "C  1": '#9467bd', "Si 1": '#8c564b', "line": '#e377c2', "O  4": '#7f7f7f', "O  3": '#bcbd22', "O  2": '#17becf', "O  1": "k", "Al 1": "b"}
replacements = {"Hlin": "H line"}

for name in cooling:
    cooling[name] = np.array(cooling[name])
    if np.max(cooling[name][:,1]) > 0.15:
        label = name
        if name in replacements:
            label = replacements[name]
        plt.semilogx(cooling[name][:,0] / Rp, cooling[name][:,1], label=label)#, color=colors[name])
        
plt.legend(frameon=False)
plt.xlabel("Radius (Rp)", fontsize=12)
plt.ylabel("Heating contribution", fontsize=12)
plt.xlim(1, 15)
plt.ylim(0, 1)
plt.figure()


#plt.loglog(radii, cooling_rate, label="Cooling rate")
plt.loglog(radii, cooling_rate / heating_rate)
#plt.legend()
plt.show()    
