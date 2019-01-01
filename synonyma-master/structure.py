
import numpy as np

hkl = np.array([
    [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1],
    [1, 0, -1], [1, -1, 0], [0, 1, -1], [2, 1, 0], [2, 0, 1], [0, 2, 1],
    [1, 2, 0], [0, 1, 2], [1, 0, 2], [1, 3, 0], [3, 1, 0], [1, 0, 3],
    [3, 0, 1], [0, 3, 1], [0, 1, 3], [3, 2, 0], [3, 0, 2], [2, 3, 0],
    [2, 0, 3], [0, 2, 3], [0, 3, 2], [4, 1, 0], [4, 0, 1], [1, 4, 0],
    [1, 0, 4], [0, 1, 4], [0, 4, 1], [5, 0, 1], [5, 1, 0], [1, 0, 5]
])

hklt = np.array([
    [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1]
])

# given coordinates, compute structure factors and intensities
def structure_factor(x):
    y = np.copy(x)
    y.shape = (-1, 3)
    sf = 1+np.sum(np.prod(y[:,np.newaxis]**hkl, axis=2), axis=0)
    return sf, sf * np.conj(sf)

class Structure:
    def __init__(self, N):
        self.N = N
        self.n = 3 * (N-1)
        s = np.random.rand(self.n)
        self.x = np.exp(2 * np.pi * 1j * s)
        self.F, self.I_F = structure_factor(self.x)

        self.random_start()

    # generate a random structure to construct a start system from
    def random_start(self):
        s = np.random.rand(self.n)
        self.x0 = np.exp(2 * np.pi * 1j * s)
        self.G, self.I_G = structure_factor(self.x0)

    # return a function that computes the homotopy
    def h_crystal(self):
        def h(x, t):
            _, I = structure_factor(x)
            return (I - (1-t) * self.I_G - t * self.I_F)[:self.n]
        return h

    # return a function that computes derivatives of homotopy
    def dh_crystal(self):
        idx = hkl[:self.n]
        def dh(i, x, t):
            y = np.copy(x)
            y.shape = (-1, 3)
            F, _ = structure_factor(x)
            d1 = idx[:,i%3] * np.prod(y[i/3]**idx, axis=1) / y[i/3,i%3]
            d2 = -idx[:,i%3] * np.prod(y[i/3]**-idx, axis=1) / y[i/3,i%3]
            return d1 * np.conj(F[:self.n]) + F[:self.n] * d2
        return dh
