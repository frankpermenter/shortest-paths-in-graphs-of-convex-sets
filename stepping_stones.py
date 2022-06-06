#!/usr/bin/env python
# coding: utf-8

# In[ ]:



# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from spp.convex_sets import Singleton, Polyhedron, CartesianProduct
from spp.convex_functions import SquaredTwoNorm
from spp.pwa_systems import PieceWiseAffineSystem, ShortestPathRegulator


# In[ ]:


# initial state
z1 = np.array([-3.5, .5, 0, 0])
q1 = z1[:2]

# target set
zK = np.array([3.5, 6.5, 0, 0])
zK = z1
qK = zK[:2]
Z = Singleton(zK)

# time horizon
K = 2

# cost matrices
q_dot_cost = .2 ** .5
Q = np.diag([1, 1, q_dot_cost, q_dot_cost])
R = np.eye(2)
S = Q # ininfluential
cost_matrices = (Q, R, S)

# In[ ]:

# configuration bounds
Dq = [
    Polyhedron.from_bounds([-4, 0], [3, 1]),
    Polyhedron.from_bounds([-6, 1], [-5, 3]),
    Polyhedron.from_bounds([4, 1], [5, 2]),
    Polyhedron.from_bounds([-4, 3], [4, 4]),
    Polyhedron.from_bounds([-5, 5], [-4, 6]),
    Polyhedron.from_bounds([5, 4], [6, 6]),
    Polyhedron.from_bounds([-3, 6], [4, 7])
]

# velocity bounds
qdot_max = np.ones(2) * 1
qdot_min = - qdot_max
Dqdot = Polyhedron.from_bounds(qdot_min, qdot_max)

# control bounds
u_max = np.ones(2) * 1
u_min = - u_max
Du = Polyhedron.from_bounds(u_min, u_max)

# pwa domains
domains = [CartesianProduct((Dqi, Dqdot, Du)) for Dqi in Dq]


# In[ ]:


# dynamics
A = np.array([
    [1, 0, 1, 0],
    [0, 1, 0, 1],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])
B = np.vstack((np.zeros((2, 2)), np.eye(2)))
Bred = B / 10
c = np.zeros(4)
dynamics = [(A, Bred, c) if i in [1, 5] else (A, B, c) for i in range(len(domains))]

# pieceiwse affine system
pwa = PieceWiseAffineSystem(dynamics, domains)


# In[ ]:


# solve optimal control problem
relaxation = 1
reg = ShortestPathRegulator(pwa, K, z1, Z, cost_matrices, relaxation=relaxation)
sol = reg.solve()
print('Cost:', sol.spp.cost)
print('Solve time:', sol.spp.time)

# unpack result
q = sol.z[:, :2]
u = sol.u

