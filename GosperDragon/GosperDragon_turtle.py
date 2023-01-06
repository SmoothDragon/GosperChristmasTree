#!/usr/bin/env python

import turtle

def gosper_curve(order, size, is_A = True):
    if(order == 0): 
        turtle.forward(size) 
        return

    for op in "A--ABA--AB++B++" if is_A else "--A--AB++BAB++B":
        gosper_op_map[op](order-1, size)
        
gosper_op_map = {"A": lambda o, size : gosper_curve(o, size, True),
                 "B": lambda o, size : gosper_curve(o, size, False),
                 "-": lambda o, size : turtle.right(60),
                 "+": lambda o, size : turtle.left(60)}

size = 3
order = 4
turtle.speed("fastest")
gosper_curve(order, size)
