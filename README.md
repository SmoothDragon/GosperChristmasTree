# GosperChristmasTree
3D Christmas tree  design based on Gosper Dragon fractal

# Idea 
Make a Christmas tree from succesive iterations of Gosper Hexagon.

# Progress
I tried to scale the shape to a point first, but the edges got so close that it could not be sliced correctly in vase mode.

Next idea is to have a different iteration at each level, and then connect the levels together.

First step is to make a routine that connects two layer that are a multiple apart in their number of points.
1 - Start with a triangle array. 
2 - Compute new array where each side replaced by _/\_
3 - Repeat for desired iterations
4 - Build a polyhedron by transitioning between each array level

* Transitions could have constant curve length which would make each level smaller, eventually into a point.
* Transistions could open up to full fractal. Maybe make each layer half as tall as the preceding layer.

# Dependencies
Listed in `requirements.txt`. Install using:
```
pip install -r requirements.txt
```
