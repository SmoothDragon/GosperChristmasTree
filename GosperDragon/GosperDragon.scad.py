#!/usr/bin/env python3

import matplotlib
# matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from solid import *
from solid.utils import *
import stl
import sys



def xy(vec):
    return vec.real, vec.imag


def gosper_iter(v0):
    r'''Apply Gosper iteration to edges.
    Sequence of complex edges ___ -> 
       _
       \
    /___\
    '''
    w = (1 + 1j*np.sqrt(3))/2  # sixth root of unity
    v1 = np.roll(v0, -1)  # All elements of v shifted to the left
    V = ((5-1j*np.sqrt(3))/(2*np.sqrt(7)))*(v1-v0)
    V = V/np.abs(V)
    M = np.vstack([
        v0,
        v0+w**4*V,
        v0+w**4*V+V,
        v0+w**4*V+2*V,
        v0+w**4*V+3*V,
        v0+2*V,
        v0+w*V+V,
        v0+w*V+2*V,
        ])
    # return entire sequence in order, each edge at a time.
    return M.flatten(order='F')  # Fortran order to read down columns

def gosper_path(p0, p1, order):
    '''
    Gosper Dragon
    A--ABA--AB++B++
    --A--AB++BAB++B
    '''
    if order == 0:
        yield p0
        # yield p1
    else:
        v = ((5-1j*np.sqrt(3))/(2*np.sqrt(7)**2))*(p1-p0)  # Lattice basis vector
        w = v*(-1 + 1j*np.sqrt(3))/2  # Second Eisenstein basis vector
        yield from gosper_path(p0, p0-v-w, order-1)
        yield from gosper_path(p0-v-w, p0-w, order-1)
        yield from list(gosper_path(p0+v-w, p0-w, order-1))[::-1]
        yield from gosper_path(p0+v-w, p0+2*v-w, order-1)
        yield from gosper_path(p0+2*v-w, p0+2*v, order-1)
        yield from list(gosper_path(p0+2*v+w, p0+2*v, order-1))[::-1]
        yield from list(gosper_path(p0+3*v+w, p0+2*v+w, order-1))[::-1]

limit = 5
# Equilateral triangle
theta = np.linspace(0, 2*np.pi, 4)[:-1]
theta = np.exp(1j*theta)
level = [np.exp(1j*theta)]
for i in range(limit-1):
    level.append(gosper_iter(level[i]))


fig, ax = plt.subplots()
# plt.axis('off')
resolution = 512

k = 5
limit = k+2
bound = (-limit, limit)
ax.set_aspect('equal')
ax.set_xlim(*bound)
ax.set_ylim(*bound)


iter = 2
dragon, = ax.plot(*xy(np.array(list(gosper_path(*theta[:2], iter)))), 'b', linewidth=1)
dragon, = ax.plot(*xy(theta[1]*np.array(list(gosper_path(*theta[:2], iter)))), 'r', linewidth=1)
dragon, = ax.plot(*xy(theta[2]*np.array(list(gosper_path(*theta[:2], iter)))), 'g', linewidth=1)
plt.show()

dragon = np.array(list(gosper_path(*theta[:2], iter)))
dragon = np.concatenate([dragon, theta[1]*dragon, theta[2]*dragon])

gosper = polygon(points=list(zip(*xy(50*dragon))))
gosper = linear_extrude(height=50)(gosper)
print(scad_render(gosper))

if __name__ == '__main__':
    import doctest
    doctest.testmod()
