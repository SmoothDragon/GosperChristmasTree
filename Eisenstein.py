#!/usr/bin/env python3

'''
Make a triangular lattice based on the Eisenstein integers.

\omega is the third root of unity
w = (-1+sqrt(3)*1j)/2


'''

import numpy as np


def xy(vec):
    return vec.real, vec.imag

def parallelogram(R, n):
    '''
    R = side length
    n = number of edges along length of R
    total of n+1 points
    '''
    w = (-1+np.sqrt(3)*1j)/2
    v = w*R/n
    L = [np.linspace(0, R, (n+1), dtype=np.complex128)]
    for i in range(n):
        L.append(L[-1]+v)
    return np.array(L)

def hexagon(R, n):
    '''Hexagon of radius R
    R = side length
    n = number of edges along length of R
    total of n+1 points

    Given a parallelogram, we rotate it three times and union to 
    get the entire hexagon.
    To avoid double counting, one of the central spines should be removed.

    Will be missing (0,0)
    '''
    w = (-1+np.sqrt(3)*1j)/2
    L = parallelogram(R, n)[:,1:]
    return np.vstack([L,w*L, w*w*L])
    # L.append(w*L[-1])
    # L.append(w*L[-1])
    # return np.array(L)



if __name__ == '__main__':
    L = (parallelogram(4,4))
    print(L)
    print()
    # print(L[:,1:])
    H = hexagon(4,10)
    print(H)

    import matplotlib.pyplot as plt
    plt.axis('equal')
    for line in H:
        plt.scatter(*xy(line))
    plt.show()
