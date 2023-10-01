import numpy as np

elements = ["H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Fe"]
weights = np.array([1, 4, 12, 14, 16, 20, 24, 28, 32, 56])
abunds = 10**np.array([12, 10.93, 8.43, 7.83, 8.69, 7.93, 7.60, 7.51, 7.12, 7.50])
abunds[2:] *= 100
print("mu", np.sum(weights * abunds) / np.sum(abunds))
#print(abunds[0] / np.sum(abunds))
print("Number fracs", abunds / np.sum(abunds))
print("Weight fracs", weights * abunds / np.sum(weights * abunds))


