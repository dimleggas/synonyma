import numpy as np

X = np.loadtxt('N=3--rec.csv', delimiter=',')
runs = np.sort(X[:,0])
phys = np.sort(X[:,1])

print runs.min()
print runs[125]
print runs[250]
print runs[375]
print runs.max()

print runs.mean()
print phys.sum() / runs.sum()
