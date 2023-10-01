import astropy.io.fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter
from scipy.interpolate import interp1d

PARSEC = 3.086e18
SEMIMAJOR_AXIS = 0.01544 * 1.5e13
DIST = 12.6 * PARSEC
CM_TO_ANGSTROM = 1e8
SPEED_OF_LIGHT = 3e10

hdul = astropy.io.fits.open("/home/stanley/packages/helium_escape/55cnc_spectrum.fits")
wavelengths = hdul[1].data["WAVELENGTH"] / CM_TO_ANGSTROM
fluxes = hdul[1].data["FLUX"] * CM_TO_ANGSTROM
print(len(fluxes))

'''interpolator = interp1d(wavelengths, fluxes)
uniform_wavelengths = np.arange(round(np.min(wavelengths)),
                                round(np.max(wavelengths)), 1)
uniform_fluxes = np.zeros(len(uniform_wavelengths))
for i, w in enumerate(uniform_wavelengths):
    cond = np.logical_and(wavelengths > w - 0.5, wavelengths < w + 0.5)
    #print("Number: ", w, np.sum(cond))
    if np.sum(cond) == 0:
        uniform_fluxes[i] = interpolator(w)
    else:
        uniform_fluxes[i] = np.sum(fluxes[cond])
        

uniform_wavelengths /= CM_TO_ANGSTROM
uniform_fluxes *= CM_TO_ANGSTROM
photon_spectrum = uniform_fluxes / (6.626e-27 * 3e10 / uniform_wavelengths)
print("Number of photons per second per cm^2:", (DIST/SEMIMAJOR_AXIS)**2 * np.sum(photon_spectrum) / CM_TO_ANGSTROM)'''

frequencies = SPEED_OF_LIGHT / wavelengths
spec_at_planet = fluxes * (DIST/SEMIMAJOR_AXIS)**2 * wavelengths**2 / SPEED_OF_LIGHT

print("interpolate")
for i in range(len(frequencies)-1, -1, -1):
    print("continue {:.6e} {:.6e}".format(
        np.log10(frequencies[i]),
        np.log10(spec_at_planet[i])))
print("luminosity solar 0.59 linear")
plt.figure()

plt.loglog(uniform_wavelengths * 1e4, uniform_fluxes)
plt.show()
