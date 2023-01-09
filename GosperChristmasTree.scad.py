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

def RecurseGosper(G):
    '''Take a gosper path to the next iteration.
    Input:  np.array(int8)
    Output: np.array(int8)
    The lower three bits of array elements contain a value modulo 6.
    ii in 0..5 represent the directions of the six roots of unity. exp(ii * 2*pi*1j/6)
    Byte stores B flag, zero bit, value mod 6: BB0321
    0x30=96 is the magic value to XOR to indicate B that won't mess up (%6)

    Gosper Dragon
    A--ABA--AB++B++
    --A--AB++BAB++B
    '''
    AB = G>>6  # 0 if A edge, 1 if B edge
    # print(AB)
    return np.stack( (  # recursion
        (G+    AB)%6,       #   A or  +A
        (G-1+  AB)%6+96,    #  -B or  -B
        (G-3+3*AB)%6+96,    # --B or   B
        (G-2     )%6+AB*96, #  +A or --B
        (G  -3*AB)%6,       # ++A or  -A
        (G  -1*AB)%6,       #   A or ++A
        (G+1-1*AB)%6+96,    #  +B or  +B
        ), 1).ravel()

def RecurseGosper2(G):
    '''Take a gosper path to the next iteration.
    Input:  np.array(int8)
    Output: np.array(int8)
    The lower three bits of array elements contain a value modulo 6.
    ii in 0..5 represent the directions of the six roots of unity. exp(ii * 2*pi*1j/6)
    Byte stores B flag, zero bit, value mod 6: BB0321
    0x30=96 is the magic value to XOR to indicate B that won't mess up (%6)

    Gosper Dragon
    A--ABA--AB++B++
    --A--AB++BAB++B
    '''
    AB = G>>6  # 0 if A edge, 1 if B edge
    # print(AB)
    return np.stack( (  # recursion
        (G+ -2*AB)%6,       #   A or --A
        (G-2-2*AB)%6,       # --A or --A
        (G-2-2*AB)%6+96,    #   B or   B
        (G-2     )%6+AB*96, #   A or ++B
        (G-4+2*AB)%6,       # --A or   A
        (G-4+2*AB)%6+96,    #   B or   B
        (G-2+2*AB)%6+96,    # ++B or ++B
        ), 1).ravel()

def show_turtle(G, size=5, last=0):
    '''Draw turtle walk on triangular lattice given a sequence of directions in 0..5
    '''
    import turtle

    turtle.speed("fastest")
    for val in list(G%6):
        turtle.left((val-last)*60)
        turtle.forward(size)
        last = val

def show_plot(G):
    fig, ax = plt.subplots()
    # plt.axis('off')
    resolution = 512

    ax.set_aspect('equal')
    points = np.cumsum(np.exp(G*2j*np.pi/6))
    ax.fill(*xy(points),'g')
    plt.show()

def show_plots(G1,G2):
    fig, axs = plt.subplots(2)
    # plt.axis('off')
    resolution = 512

    axs[0].set_aspect('equal')
    axs[1].set_aspect('equal')
    points = np.cumsum(np.exp(G1*2j*np.pi/6))
    axs[0].fill(*xy(points),'g')
    points = np.cumsum(np.exp(G2*2j*np.pi/6))
    axs[1].fill(*xy(points),'b')
    plt.show()


def GosperTrim(G):
    result = [G[0]]
    for ii in range(1, len(G)):
        if (G[ii-1]-G[ii])%6 == 4:
            result[-1] = (G[ii]-1)%6
        else:
            result.append(G[ii])
    if (G[-1]-G[0])%6 == 4:
        result[0] = (G[0]-1)%6
        result.pop()
        # result = result[:-1]
    print(G)
    result = np.array(result, dtype=np.int8)
    print(result)
    return result



v = (np.array([0,2,4], dtype=np.int8))
# v = (np.array([0,1+96,2,3+96,4, 5+96], dtype=np.int8))
# v = RecurseGosper(np.array(range(6)))
# v = (np.array([4,2,0]))
# v = RecurseGosper2(v)
# v = RecurseGosper2(v)
# v = RecurseGosper2(v)
v = RecurseGosper(v)
v = RecurseGosper(v)
v = RecurseGosper(v)
v = RecurseGosper(v)

print(v)
v=v%6
# show_turtle(v); input()
# show_plot(v)
w = GosperTrim(v)
n=1
while len(w):
    show_plots(v,w)
    v,w = w,GosperTrim(w)
    n += 1
print(f'Iterations: {n}')


exit(0)

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
limit = 2
bound = (-limit, limit)
ax.set_aspect('equal')
ax.set_xlim(*bound)
ax.set_ylim(*bound)


iter = 3
# dragon, = ax.plot(*xy(np.array(list(gosper_path(*theta[:2], iter)))), 'b', linewidth=1)

# plt.show()

dragon0 = np.array(list(gosper_path(*theta[:2], 0)))
dragon1 = np.array(list(gosper_path(*theta[:2], 1)))
dragon = np.array(list(gosper_path(*theta[:2], iter)))
dragon = np.concatenate([dragon, theta[1]*dragon, theta[2]*dragon])
# dragon0 = np.concatenate([dragon0, theta[1]*dragon0, theta[2]*dragon0])

ax.fill(*xy(dragon),'g')
# ax.plot(*xy(dragon),'r')
# ax.plot(*xy(dragon0),'g')
# ax.plot(*xy(dragon1),'b')
# print(dragon0)
# plt.show()
# exit(0)

diameter = 100
height = 100
edge = diameter/5
stud_shift = edge/2
stud_dia = 2.5
limit = 5


gosper = polygon(points=list(zip(*xy(edge*dragon))))
# gosper = linear_extrude(height=2*height, center=True, slices=10, twist=100)(gosper)
gosper = linear_extrude(height=2*height, center=True)(gosper)
# stud_hole = cylinder(d=stud_dia, h=2*height, center=True, segments=8)
# stud_holes = union()([translate([stud_shift*(i+(j&1)/2), stud_shift*j*np.sqrt(3)/2,0])(stud_hole) for i in range(-limit, limit+1) for j in range(-limit, limit+1)])
# stud_holes = translate([-edge/4,edge*np.sqrt(3)/12,0])(stud_holes)
# stud_holes = rotate([0,0,-8])(stud_holes)
# stud_holes = intersection()(stud_holes,
                            # cylinder(d=.97*diameter, h=4*height, center=True, segments=6))
# gosper = rotate([0,0,8])(gosper)
base = cylinder(d1=diameter, d2=0, h=height, center=True, segments=256)
gosper = intersection()(gosper,base)
print(scad_render(gosper))

if __name__ == '__main__':
    import doctest
    doctest.testmod()
