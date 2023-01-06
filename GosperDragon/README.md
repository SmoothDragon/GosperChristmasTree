Install:
pip install solid solid-toolbox

Idea is to make a Christmas tree design from a space filling curve.
Space filling curves don't tend to be loops.
I liked the Gosper Curve, but alas, not a loop.
I discoverd that three Gosper dragons can be joined to make a loop which is a Gosper curve.
So I should be able to make a closed loop for openscad, extrude it upward and then intersect the loop with a cone to make a Christmas tree like shape that still can be printed in vase mode.

Turtle sequence for Gosper Curve

A-B--B+A++AA+B-
+A-BB--B-A++A+B

Gosper Dragon
A--ABA--AB++B++
--A--AB++BAB++B

Triple flowsnake from:
https://blogs.scientificamerican.com/roots-of-unity/a-few-of-my-favorite-spaces-space-filling-curves/
