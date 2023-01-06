#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
plt.axis('off')
resolution = 256

num = 5
nloops = 1  # nloops*num should be an integer
r = 1/num
k = 5

ax.set_aspect('equal')
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)

# big circle:
theta = np.linspace(0, 2*np.pi*nloops, resolution*nloops+1)
x = np.cos(theta)
y = np.sin(theta)

# hypocycloid
hax = (k-1)*np.cos(theta) + np.cos((k-1)*theta)
hay = (k-1)*np.sin(theta) - np.sin((k-1)*theta)

hbx = k*np.cos(theta) + np.cos(k*theta)
hby = k*np.sin(theta) - np.sin(k*theta)

# Rolling circle
sx = r*np.cos(theta)
sy = r*np.sin(theta)

big_circle, = ax.plot(x, y, 'b', linewidth=2)
dot, = ax.plot([], [], 'ko', ms=10)
hypocycloid_a, = ax.plot([], [], 'r-', linewidth=3.3)
hypocycloid_b, = ax.plot([], [], 'g-', linewidth=3.3)
small_circle, = ax.plot([], [], 'k-', linewidth=2)
line, = ax.plot([], [], 'k-', linewidth=2)


def frame_generator(step=8):
    # Drawing path
    for i in range(0, len(theta), step):
        X, Y = x[i], y[i]  # Center
        dot.set_data([X], [Y])
        h3x = np.cos(-theta[i]/k)*(hax) - np.sin(-theta[i]/k)*(hay) + X
        h3y = np.sin(-theta[i]/k)*(hax) + np.cos(-theta[i]/k)*(hay) + Y
        # hypocycloid_a.set_data(hax+X, hay+Y)
        hypocycloid_a.set_data(h3x, h3y)
        hypocycloid_b.set_data(hbx, hby)
        yield


def no_op(*args, **kwargs):
    return


ani = FuncAnimation(fig, no_op, frames=frame_generator(), save_count=300, blit=False)

ani.save('hypo3to4.gif', writer='imagemagick', fps=10, dpi=75)
# plt.show()

"""
class Hypocycloid:

    def __init__(self, ratio=3, frames=100, ncycles=1):
        self.frames = frames
        self.ncycles = ncycles
        self.fig, self.ax = plt.subplots()
        self.ax.set_aspect('equal')
        plt.axis('off')


        ##big circle:
        theta = np.linspace(0, 2*np.pi, 100)
        x = np.cos(theta)
        y = np.sin(theta)

        # self.white_circle, = self.ax.plot(1.0001*x, 1.0001*y, 'w-', linewidth=4)
        self.big_circle, = self.ax.plot(x, y, 'b', linewidth=2)
        # self.dot_circle, = self.ax.plot(x*(1-1/ratio), y*(1-1/ratio), color='gray', linestyle='dashed', linewidth=2)
        self.dot_circle, = self.ax.plot(x*(1-1/ratio), y*(1-1/ratio), 'k--', linewidth=1.5)

        ##small circle:
        self.small_r = 1./ratio
        r = self.small_r
        x = r*np.cos(theta)+1-r
        y = r*np.sin(theta)
        self.small_circle, = self.ax.plot(x,y,'k-',linewidth=2)
        self.small_circle2, = self.ax.plot(-x,-y,'k-',linewidth=2)

        ##line and dot:
        self.line, = self.ax.plot([1-r,1],[0,0],'k-',linewidth=2)
        self.line2, = self.ax.plot([],[],'k-',linewidth=2)
        self.dot, = self.ax.plot([1-r],[0], 'ko', ms=10)
        self.dot2, = self.ax.plot([1-r],[0], 'ko', ms=10)
        ##hypocycloid:
        self.hypocycloid, = self.ax.plot([],[],'r-',linewidth=3.3)
        self.hypocycloid2, = self.ax.plot([],[],'r-',linewidth=3.3)


        self.animation = FuncAnimation(
            self.fig, self.animate,
            frames=self.frames*self.ncycles,
            interval=50, blit=False,
            repeat_delay=2000,
        )

    def update_small_circle(self, phi):
        theta = np.linspace(0,2*np.pi,100)
        x = self.small_r*np.cos(theta)+(1-self.small_r)*np.cos(phi)
        y = self.small_r*np.sin(theta)+(1-self.small_r)*np.sin(phi)
        self.small_circle.set_data(x,y)
        self.small_circle2.set_data(-x,-y)


    def update_hypocycloid(self, phis):
        r = self.small_r
        x = (1-r)*np.cos(phis)+r*np.cos((1-r)/r*phis)
        y = (1-r)*np.sin(phis)-r*np.sin((1-r)/r*phis)
        self.hypocycloid.set_data(x,y)
        self.hypocycloid2.set_data(y,x)

        center = [(1-r)*np.cos(phis[-1]), (1-r)*np.sin(phis[-1])]

        self.line.set_data([center[0],x[-1]],[center[1],y[-1]])
        self.line2.set_data([center[0],-x[-1]],[center[1],-y[-1]])
        self.dot.set_data([center[0]], [center[1]])
        self.dot2.set_data([-center[0]], [-center[1]])

    def animate(self, frame):
        frame = frame+1
        phi = 2*np.pi*frame/self.frames
        self.update_small_circle(phi)
        self.update_hypocycloid(np.linspace(0,phi,frame))

hypo = Hypocycloid(ratio=5, frames=40, ncycles=1)

##un-comment the next line, if you want to save the animation as gif:
hypo.animation.save('hypo.gif', writer='imagemagick', fps=10, dpi=75)

# plt.show()
"""
