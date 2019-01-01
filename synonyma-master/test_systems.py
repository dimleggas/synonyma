
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
