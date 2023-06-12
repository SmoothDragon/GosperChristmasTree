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
import solid as sd

import tempfile
import subprocess


def scad2svg(obj):
    with tempfile.NamedTemporaryFile(mode='w', dir='.') as tmp:
        tmp.writelines(sd.scad_render(obj))
        tmp.flush()
        command = f'openscad --export-format svg -o - {tmp.name}'
        p = subprocess.run(command, capture_output=True, shell=True, encoding='utf-8')
    return p.stdout

def setSVGproperty(svg, stroke=None, fill=None, stroke_width=None):
    if fill is not None:
        svg = svg.replace('fill="lightgray"', f'fill="{fill}"')
    if stroke is not None:
        svg = svg.replace('stroke="black"', f'stroke="{stroke}"')
    if stroke_width is not None:
        svg = svg.replace('stroke-width="0.5"', f'stroke-width="{stroke_width}"')
    return svg

def joinSCAD_SVG(svg0,svg1):
    result = svg0.split('\n')[:-2]+svg1.split('\n')[4:]
    return '\n'.join(result)

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
# dragon, = ax.plot(*xy(np.array(list(gosper_path(*theta[:2], iter)))), 'b', linewidth=1)
# dragon, = ax.plot(*xy(theta[1]*np.array(list(gosper_path(*theta[:2], iter)))), 'r', linewidth=1)
# dragon, = ax.plot(*xy(theta[2]*np.array(list(gosper_path(*theta[:2], iter)))), 'g', linewidth=1)
# plt.show()

dragon = np.array(list(gosper_path(*theta[:2], iter)))
dragon = np.concatenate([dragon, theta[1]*dragon, theta[2]*dragon])

diameter = 100
height = 3
edge = diameter/5
stud_shift = edge/2
stud_dia = 2.5
limit = 5


gosper = polygon(points=list(zip(*xy(edge*dragon))))
r = 4
gosper = offset(r)(gosper)
gosper = offset(-2*r)(gosper)
gosper = offset(r)(gosper)
stud_hole = circle(d=stud_dia,  segments=8)
stud_holes = union()([translate([stud_shift*(i+(j&1)/2), stud_shift*j*np.sqrt(3)/2])(stud_hole) for i in range(-limit, limit+1) for j in range(-limit, limit+1)])
stud_holes = translate([-edge/4,edge*np.sqrt(3)/12])(stud_holes)
# stud_holes = rotate([0,0,-8])(stud_holes)
stud_holes = intersection()(stud_holes,
                            circle(d=.97*diameter, segments=6))
# gosper = rotate([0,0,8])(gosper)
gosper = rotate(8)(gosper)
base = circle(d=.6*diameter,  segments=6)
base = offset(.2*diameter)(base)
# base -= gosper
# base -= stud_holes
# print(scad_render(base-stud_holes))
svg_bb = scad2svg(base)
svg_bb = setSVGproperty(svg_bb, fill='white')
svg_gosper = scad2svg(gosper)
svg_gosper = setSVGproperty(svg_gosper, fill='red')
svg_holes = scad2svg(stud_holes)
svg_holes = setSVGproperty(svg_holes, fill='black')
total = joinSCAD_SVG(svg_bb, svg_gosper)
total = joinSCAD_SVG(total, svg_holes)
print(total)

# if __name__ == '__main__':
    # import doctest
    # doctest.testmod()

    # svg_back = scad2svg(backdrop(x,y,square))
    # svg_back = setSVGproperty(svg_back, fill='white')
    # svg_board = scad2svg(board(x,y,square))
    # svg_board = setSVGproperty(svg_board, fill='black', stroke_width=0)
    # svg_back = joinSCAD_SVG(svg_bb, svg_back)
    # total = joinSCAD_SVG(svg_back, svg_board)
    # print(total)

