import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_iters(depths, quantity, starts, title=""):
    plt.figure()
    for i in range(len(starts) - 1):
        #if i != len(starts) - 2: continue
        if i % 10 != 0: continue
        first = starts[i]
        end = starts[i+1]
        plt.semilogy(radii[first:end], 10**quantity[first:end])
    plt.title(title)


#Rp = 1.3e9
R_earth = 6.378e8
Rp = 2.15 * R_earth

data = np.loadtxt(sys.argv[1])
#data = data[-int(len(data)/2):]
#data = data[0:492]
#for iter in range(30):
#    print(data[4


depths, Te, Htot, hden, eden, H_molec_ratio, HI, HII, HeI, HeII, HeIII = data.T[0:11]

'''if np.max(depths) / Rp + 1 <= 12:
    max_R = 11.92451
elif np.max(depths) / Rp + 1 > 14:
    max_R = 14.88
else:
    assert(False)'''

radii = 14.9 - depths/Rp #not exact, found by trial and error
#print(depths)
starts = np.where(np.abs(depths - 0.25) < 1e-3)[0]
starts = np.append(starts, len(depths))
print(starts)
plot_iters(radii, Te, starts, "Te")
plot_iters(radii, HI, starts, "HI")
plot_iters(radii, HII, starts, "HII")
plot_iters(radii, HeI, starts, "HeI")
plot_iters(radii, HeII, starts, "HeII")
plot_iters(radii, HeIII, starts, "HeIII")

plot_iters(radii, Htot, starts, "Htot")
plot_iters(radii, hden, starts, "hden")
plot_iters(radii, eden, starts, "eden")

plt.figure()
last = 900
plt.plot(radii[-last:], HI[-last:])
plt.plot(radii[-last:], HII[-last:])
plt.plot(radii[-last:], HeI[-last:])
plt.plot(radii[-last:], HeII[-last:])
plt.plot(radii[-last:], eden[-last:])
plt.show()
