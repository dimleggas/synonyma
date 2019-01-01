
import numpy as np
from scipy.spatial.distance import pdist, squareform

'''
    structures; homotopy and its jacobian
'''
hkl = np.array([
    [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 0, -1], [1, -1, 0], 
    [0, 1, -1], [2, 1, 0], [2, 0, 1], [0, 2, 1], [1, 2, 0], [0, 1, 2], [1, 0, 2], [1, 3, 0], 
    [3, 1, 0], [1, 0, 3], [3, 0, 1], [0, 3, 1], [0, 1, 3], [3, 2, 0], [3, 0, 2], [2, 3, 0],
    [2, 0, 3], [0, 2, 3], [0, 3, 2], [4, 1, 0], [4, 0, 1], [1, 4, 0],
    [1, 0, 4], [0, 1, 4], [0, 4, 1], [5, 0, 1], [5, 1, 0], [1, 0, 5]
    ])

def structure(N):
    return np.exp(2 * np.pi * 1j * np.random.rand(N-1, 3))

def real_structure(N):
    d = 1. / np.cbrt(N)
    delta = 0.1 * d
    real = False
    while not real:
        s = np.random.rand(N-1, 3)
        x = np.r_[np.zeros((1, 3)), s]
        D = squareform(pdist(x))
        real = True
        for row in D:
            if np.any(row[row>0] < d-delta) or np.all(row[row>0] > d+delta):
                real = False
    return np.exp(2 * np.pi * 1j * s)

def structure_factor(N, x):
    p = hkl[:3*(N-1)].reshape(3*(N-1), 1, 3)
    F1 = 1 + (x**p).prod(-1).sum(-1)
    F2 = 1 + (x**(-p)).prod(-1).sum(-1)
    return F1, F2, F1*F2

def homotopy_stuff(I0, I1, N):
    n = 3*(N-1)
    def _homotopy(x, t):
        F1, F2, I = structure_factor(N, x)
        return I - (1-t)*I0 - t*I1, F1, F2
    def _jacobian(x, F1, F2):
        n = 3*(N-1)
        J = np.zeros((n,n), dtype=np.complex64)
        for i in range(3*(N-1)):
            for j in range(3*(N-1)):
                p1, p2 = np.copy(hkl[i]), -np.copy(hkl[i])
                c1, c2 = p1[j%3], p2[j%3]
                p1[j%3] -= 1
                p2[j%3] -= 1
                x_ = np.copy(x[j/3])
                dF1, dF2 = c1*(x_**p1).prod(), c2*(x_**p2).prod()
                J[i, j] = dF1 * F2[i] + F1[i] * dF2  
        return np.linalg.inv(J)
    return _homotopy, _jacobian

'''
    homotopy path
'''
def newton(x, t, homotopy, jacobian):
    N = x.shape[0] + 1
    tol, step, maxit, it = 1e-10, 100, 50, 0
    while it < maxit and tol < step:
        it += 1
        H, F1, F2 = homotopy(x, t)
        J = jacobian(x, F1, F2)
        x = x - np.dot(J, H).reshape(N-1, 3)
        step = np.sqrt(np.sum(np.square(H - homotopy(x, t)[0])))
    return x

def homotopy(N, x0, I1, dt=.05):
    _, _, I0 = structure_factor(N, x0)
    homo, jaco = homotopy_stuff(I0, I1, N)
    x = np.copy(x0)
    for t in np.arange(dt, 1+dt, dt):
        x = newton(x, t, homo, jaco)
    return x

def sort(x):
    return x[np.argsort(x[:,0])]

def same_structure_Q(x, x_target):
    N = x.shape[0] + 1
    x_target = sort(x_target)
    if np.allclose(sort(x), x_target):
        return True
    if np.allclose(sort(np.conj(x)), x_target):
        return True
    for i in range(N-1):
        x_shift = np.copy(x)
        shift = np.conj(x[i])
        x_shift = x_shift * shift
        x_shift[i] = shift
        if np.allclose(sort(x_shift), x_target):
            return True
        if np.allclose(sort(np.conj(x_shift)), x_target):
            return True
    return False

def solve_structure(N, x1):
    _, _, I1 = structure_factor(N, x1)
    solved = False
    num_runs, num_phys = 0, 0
    while not solved:
        num_runs+=1
        x0 = structure(N)
        x_ = homotopy(N, x0, I1)
        if np.allclose(np.abs(x_), 1):
            num_phys+=1
            if same_structure_Q(x_, x1):
                solved = True
    return x_, num_runs, num_phys

def recursive_solve(N, I1):
    if N == 2:
        while True:
            x0 = structure(N)
            x_ = homotopy(N, x0, I1)
            if np.allclose(np.abs(x_), 1):
                return x_
    else:
        I1_ = (N-1)**2. / N**2. * I1[:-3]
        solved = False
        num_runs, num_phys = 0, 0
        while not solved:
            num_runs+=1
            x0_ = recursive_solve(N-1, I1_)
            x0 = np.r_[x0_, np.exp(2 * np.pi * 1j * np.random.rand(1,3))]
            x_ = homotopy(N, x0, I1)
            if np.allclose(np.abs(x_), 1):
                num_phys+=1
                if same_structure_Q(x_, x1):
                    solved = True
        return x_, num_runs, num_phys

N=3
c1, c2 = [], []
for i in range(500):
    x1 = structure(N)
    _, _, I1 = structure_factor(N, x1)
    _, num_runs, num_phys = recursive_solve(N, I1)
    c1.append(num_runs)
    c2.append(num_phys)
np.savetxt('N=3--rec.csv', np.c_[c1, c2].astype(np.int32), delimiter=',')
