#!/usr/bin/env python3

'''Make a Christmas Tree based on Gosper Dragon Fractal

'''


import sys
import numpy as np
import stl
import matplotlib.pyplot as plt

def parseInput():
    parser = argparse.ArgumentParser(
        description=HELP,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
    parser.add_argument(
        '-v', '--verbose',
        type=int,
        dest='verbose',
        nargs='?',
        default=logging.ERROR//10,  # Default verbosity
        help='Provide verbose output (1=DEBUG, 2=INFO, 3=WARN, 4=ERROR, 5=CRITICAL)',
        )
    parser.add_argument(
        '-V', '--version',
        action='version',
        version='%(prog)s 1.0',
        )
    parser.add_argument(
        '-n', '--number',
        type=str,
        dest='number',
        default='100',
        help='Number to convert',
        )
    parser.add_argument(
        '-p', '--precision',
        type=float,
        dest='prec',
        nargs='?',
        default=.0001,
        help='New representation will be within given precision.',
        )
    parser.add_argument(
        'args',
        nargs=argparse.REMAINDER,
        help='List of numbers to convert',
        )
    return parser.parse_args()

def setLogging(verbose):
    if verbose is None:
        level = logging.INFO  # Default verbosity level if flag is used without level
    else:
        level = verbose * 10  # Use verbosity number specified
    logging.basicConfig(level=level, format='%(levelname)-8s %(message)s')
    logging.debug('Using verbosity level {level:d}'.format(**locals()))

def plot_SVG(filename, x, y, closed=True, fill=True):
    if closed:
        x = np.concatenate((x, x[0:1]))
        y = np.concatenate((y, y[0:1]))
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    plt.axis('off')
    if fill:
        ax.fill(x, y, 'g')
    else:
        ax.plot(x, y, 'b', linewidth=2)
    plt.savefig(filename)

def xy(vec):
    return vec.real, vec.imag

def midpoint(v, alpha=.5, shift=0):
    '''Returns a vector of midpoints
    w[i] = (v[i] + v[i+1])/2
    When alpha==0, the vectors is just v.
    When alpha==1, the vectors is just v shifted.
    >>> midpoint(np.array([0,2,4]))
    array([1., 3., 2.])
    '''
    return (1-alpha)*v + alpha*np.roll(v, -1) 

def extender(v, scale=2):
    '''Increase edge list size
    The edge list increases by the scale factor by subdividing each
    edge into "scale" many edges.
    >>> extender(np.array([0,1,1j]))
    array([0. +0.j , 0.5+0.j , 1. +0.j , 0.5+0.5j, 0. +1.j , 0. +0.5j])
    '''
    w = np.roll(v, -1)  # All elements of v shifted to the left
    M = np.vstack([((scale-i)*v+i*w)/scale for i in range(scale)])
    # return entire sequence in order, each edge at a time.
    return M.flatten(order='F')  # Fortran order to read down columns

def binder_single(v,w,z0,z1):
    '''Create all the triangles joining two ringed levels.
    Let len(v)=k*len(w)
    Then k edges of v will attach to a single point of w.
    A single final triangle will connect the edge of w to the next point of v.
    '''
    m = len(w)
    n = len(v)
    k = n // m
    assert n == m*k
    # v and w are stacks of (x,y,z) points
    v0 = np.vstack([*xy(v), z0*np.ones(len(v))]).transpose()
    v1 = np.roll(v0, -1, axis=0)
    w0 = np.vstack([*xy(w), z1*np.ones(len(w))]).transpose()
    w1 = np.roll(w0, -1, axis=0)
    # Each 3xnx3 array starts as a 3 long array of stacks of (x,y,z)
    # All faces are counterclockwise
    faces = []
    for i in range(k):
        faces.append(np.array([v0[i::k], v1[i::k], w0]))
    faces.append(np.array([w0, v1[k-1::k], w1]))
    # After transpose, it becomes a nx3x3 array, or stack of triangles
    faces = [F.transpose(1,0,2) for F in faces]
    # Combine all faces into one long array of 3x3 triangle faces
    return np.vstack(faces)


def koch_iter(v):
    '''Apply Koch snowflake iteration to edges.
    Sequence of complex edges ___ -> _/\_
    Bulge sticks our to the right, so counterclockwise recommended.
    >>> koch_iter(np.array([0,1]))
    array([0.        +0.j        , 0.33333333+0.j        ,
           0.5       -0.28867513j, 0.66666667+0.j        ,
           1.        +0.j        , 0.66666667+0.j        ,
           0.5       +0.28867513j, 0.33333333+0.j        ])
    '''
    w = np.roll(v, -1)  # All elements of v shifted to the left
    # M[i] = vector of ith points (0<=i<4) in _/\_ -> 0_1/2\3_
    M = np.vstack([
        v,
        2/3*v+1/3*w,  # 1/3 of way to next point
        (v+w)/2 + 1j*(v-w)*np.sqrt(3)/6,  # point sticking out
        1/3*v+2/3*w,  # 2/3 of way to next point
        ])
    # return entire sequence in order, each edge at a time.
    return M.flatten(order='F')  # Fortran order to read down columns


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
        # v0+w*V+V,
        # v0+w*V+2*V,
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



def concave_triangulation(facet, epsilon=1e-10):
    '''Given a facet, find points that bend counter-clockwise, followed by
    points that bend clockwise. These points can be triangulated away using
    the preceding point in the path.

    This function uses the compress function in numpy to keep things speedy.
    >>> concave_triangulation(np.array([-.5,1j,-1,-1j]))
    {'faces': array([[-0. -1.j, -0.5+0.j, -1. +0.j]]), 'path': array([-0.5+0.j,  0. +1.j, -1. +0.j])}
    '''
    center = facet
    ahead = np.roll(facet, -1)
    behind = np.roll(facet, 1)
    a = ahead - center
    b = center - behind
    pos = (a*np.conj(b)).imag
    # postive bend followed by negative
    # marked = np.roll(pos<epsilon, -1) & (pos>-epsilon)
    marked = (np.roll(pos, -1)<epsilon) & (pos>epsilon)
    faces = np.array([center.compress(marked),
        ahead.compress(marked),
        behind.compress(marked),
              ]).transpose()
    path = facet.compress(~marked)
    return {'faces':faces, 'path':path}

def path_triangulation(path):
    faces = []
    count = 0
    while True:
        if is_convex(path):
            faces.append(convex_triangulation(path))
            break
        result = concave_triangulation(path)
        path = result['path']
        faces.append(result['faces'])
    # extra = np.array([path, .1*np.ones(len(path)), -.1*np.ones(len(path))]).transpose()
    # faces.append(extra)
    return np.vstack(faces)


def is_convex(facet, epsilon=1e-10):
    '''Tests facet for convexity
    Eplison should be nonzero to avoid convexity errors when points are 
    colinear. In a near colinear case, the convex triangulation should
    work anyway, but this could potentially be a problem.
    >>> is_convex(np.array([1,1j,-1,-1j]))
    True
    >>> is_convex(np.array([-.5,1j,-1,-1j]))
    False
    >>> is_convex(np.array([0,1j,-1,-1j]))
    True
    '''
    center = facet
    ahead = np.roll(facet, -1)
    behind = np.roll(facet, 1)
    a = ahead - center
    b = center - behind
    pos = (a*np.conj(b)).imag
    # print(pos, file=sys.stderr)
    return all(pos>-epsilon)

def convex_triangulation(facet):
    '''Given a convex planar facet, make triangles from the perimeter to the center
    Designed to work with complex points
    '''
    try: 
        assert(is_convex(facet))
    except:
        print(facet)
    n = len(facet)
    barycenter = sum(facet)/n
    faces = np.vstack([facet, np.roll(facet,-1), np.array([barycenter]*n)])
    return faces.transpose(1,0)

def complex_to_3d(array, height=0):
    '''Given an array of complex points, convert then to a similar array
    with the x and y points separated with z as height.
    a+bi -> (a,b,height) for all points
    '''
    M = np.array([*xy(array),height*np.ones(shape=array.shape)])
    M = M.transpose(1,2,0)
    return M

def regular_polygon(n):
    '''
    >>> np.testing.assert_almost_equal(regular_polygon(4),\
        np.array([1,1j,-1,-1j]))
    '''
    theta = np.linspace(0, 2*np.pi, n+1)[:-1]
    return np.exp(1j*theta)

def faces2mesh(faces):
    mesh = stl.Mesh(np.zeros(faces.shape[0], dtype=stl.Mesh.dtype))
    mesh.vectors = faces
    return mesh

if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # Save to outfile dropping the final ".py"
    filename = '.'.join(sys.argv[0].split('.')[:-1])

    # Equilateral triangle
    theta = np.linspace(0, 2*np.pi, 4)[:-1]
    theta = np.exp(1j*theta)
    ring = regular_polygon(3)
    # dragon = np.array(list(gosper_path(*theta[:2], iter)))
    # dragon = np.concatenate([dragon, theta[1]*dragon, theta[2]*dragon])

    ring2 = gosper_iter(ring)

    iter = 0
    fig, ax = plt.subplots()
    # dragon, = ax.plot(*xy(np.array(list(gosper_path(*theta[:2], iter)))), 'b', linewidth=1)
    # dragon, = ax.plot(*xy(theta[1]*np.array(list(gosper_path(*theta[:2], iter)))), 'r', linewidth=1)
    # dragon, = ax.plot(*xy(theta[2]*np.array(list(gosper_path(*theta[:2], iter)))), 'g', linewidth=1)
    points = np.array(list(gosper_path(*theta[:2], iter)))
    dragon = np.concatenate((points, theta[1]*points, theta[2]*points))
    plt.axis('equal')
    ax.fill(*xy(dragon))
    plt.show()

    # plot_SVG('gosper.svg', *xy(np.concatenate((ring, ring[0:1]))))
    plot_SVG('gosper.svg', *xy(ring))
    plot_SVG('gosper2.svg', *xy(ring2))
    plot_SVG('dragon.svg', *xy(dragon))
    ring = gosper_iter(ring2)
    ring, ring2 = ring2, ring
    print(len(ring))
    print(len(ring2))
    faces = binder_single(ring2, ring, 0, 10)
    # faces = binder_single(10*regular_polygon(6), 10*regular_polygon(3),0,10)
    mesh = faces2mesh(faces)
    mesh.save(filename, mode=stl.Mode.ASCII)  # For debugging
    # mesh.save(filename)

