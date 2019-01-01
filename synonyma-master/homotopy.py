
import numpy as np

class Homotopy:
	"""
	Homotopy: Single Track Polynomial Homotopy Continuation

	Numerical homotopy for calculating the solutions of an nxm system
	of polynomials equations.

	Parameters
	----------
	n : int
		number of equations
	h : function of ndarray, float returning ndarray
	 	The homotopy H(t, x; F, G) such that H(1, x)=F(x) the target system and
		H(0, x)=G(x) the start system with the same dimension as F and trivial
		solutions. e.g. the linear homotopy H(t, x; F, G) := (1-t)G(x) + tF(x)
	dh : function of int, ndarray, float returning ndarray
		The derivative of H w.r.t. the ith variable (numberical by default)
	dt : float
		timestep (0.01 by default)

	Attributes
	----------
	t : float
		time parameter in [0, 1]
	ht : ndarray of length n
		the H(x, t) evaluated at some x and t
	jt : ndarray of size nxn
		the inverse of the jacobian of H(x, t)
	x : ndarray
		solution to H(t, x) == 0

	"""

	def __init__(self, n, h, dh=None, dt=.01):
		self.n = n
		self.h = h
		self.dh = dh
		self.dt = dt

		self.t = 0
		self.ht = np.zeros(n)
		self.jt = np.zeros((n, n))
		self.x = np.zeros(n)

	# set a new homotopy definition (start system changed)
	def init_h(self, h, dh=None):
		self.h = h
		self.dh = dh

	# sets the solution to G as the solution of H
	def init_x(self, x):
		self.t = 0
		self.x = x

	# derivative estimation
	def deriv(self, i, x1, t):
		h = 1.e-10
		x2 = np.copy(x1)
		x2[i] += h
		return (self.h(x2, t) - self.h(x1, t)) / h

	# compute H at current x and t
	def homotopy(self):
		self.ht = self.h(self.x, self.t)

	# compute Jacobian of H at current x and t
	def jacobian(self):
		j = np.zeros((self.n, self.n)).astype(np.complex128)
		for i in range(self.n):
			if self.dh == None:
				j[:, i] = self.deriv(i, self.x, self.t)
			else:
				j[:, i] = self.dh(i, self.x, self.t)
		self.jt = np.linalg.inv(j)

	# the Newton corrector
	def newton(self):
		tol, maxit, it, converge = 1.e-8, 20, 0, False
		while not converge and it < maxit:
			it = it + 1
			xlast = self.x
			self.homotopy()
			self.jacobian()
			self.x = xlast - np.dot(self.jt, self.ht)
			converge = np.linalg.norm(self.x - xlast) < tol
		print it

	# Find solution of F(x)=H(x,1)==0
	def track(self):
		while self.t < 1:
			print self.t
			self.t = np.min([1, self.t + self.dt])
			self.newton()
