"""
A simulation written using the QInfo library to implement Deutsch's algorithm. 

A function f maps {0,1} to {0,1}
The function can be one of the four forms:
f(x) = 0, f(x) = 1  constant functions.
f(x) = x, f(x) = x' balanced functions.

The objective is to determine whether the function 
is a constant function or a balanced function when kept in a blackbox.

Classically we need 2 queries to determine this. 
Using quantum we need just one. 
"""

__author__ = "Madhav Jivrajani"

import numpy as np 
from QInfo import Register 
# from mapping import mapping_2 as m2
from QGates import *
import math 
import time 

I = np.array([[1,0],[0,1]])

def quantumOracle():
    """
    Returns a dictionary of different quantum oracles for the following four cases:
    const_0 :   f(x) = 0
    const_1 :   f(x) = 1
    bal_0   :   f(x) = x
    bal_1   :   f(x) = x'

    Implementing for f(x) = 1. Feel free to try out the other cases. 
    """
    const_0 = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    const_1 = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
    bal_0   = cnot
    bal_1   = np.array([[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]])

    oracles = {"const_0":const_0,"const_1":const_1,"bal_0":bal_0,"bal_1":bal_1}

    return oracles

reg = Register(2).initialize()

reg[1].pauliX()

reg[0].hadamard()
reg[1].hadamard()

#oracle
reg[1].cNOT(reg[0])

reg[0].hadamard()
reg[1].hadamard()

#Measurement
reg[0].measureNum(1)

#print out probabilities in basis states.
print(reg[0].measureRes)


