
import sys
import numpy as np
from homotopy import Homotopy
from structure import Structure

def h(x, t):
    f1 = x[1]**3-6*x[0]**2-2*x[0]*x[1]+5*x[1]**2+2*x[0]
    f2 = x[0]**2+x[1]**2-2
    g1 = x[0]**3-1
    g2 = x[1]**2-1
    return np.array([(1-t)*g1+t*f1,(1-t)*g2+t*f2])

def dh(i, x, t):
    if i == 0:
        return np.array([(1-t)*3*x[0]**2+t*(-12*x[0]-2*x[1]+2), t*(2*x[0])])
    else:
        return np.array([t*(3*x[1]**2-2*x[0]+10*x[1]),(1-t)*2*x[1]+t*(2*x[1])])

def main():

    # crys = Structure(N)
    # h = crys.h_crystal()
    # dh = crys.dh_crystal()
    homo = Homotopy(2, h, dh)
    homo.init_x([-.5 + np.divide(np.sqrt(3),2) * 1j, -1])
    print homo.x
    homo.track()
    # count = 1
    # while not np.all(np.abs(homo.x) - 1 < 1.e-4):
    #     count = count + 1
    #     crys.random_start()
    #     homo.init_h(crys.h_crystal(), crys.dh_crystal())
    #     homo.init_x(crys.x0)
    #     homo.track()
    #     print count,homo.x, np.abs(homo.x), homo.h(homo.x, 1)
    #
    # print count
    print homo.x
    # print np.abs(homo.x)
    print homo.h(homo.x, 1)

if __name__ == '__main__':
    main()
