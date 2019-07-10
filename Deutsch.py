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
import QInfo as qi 
from QGates import cnot, H
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

def quantumDeutsch():
    oracle = quantumOracle()["const_1"]

    #initialise a quantum register with 2 qubits in states |0> and |1>
    register = [qi.construct_standard_basis(1)[i] for i in range(0,2)] # The state now is |00> 

    #Creating superposition by applying hadamard transform to both qubits. 
    register[0] = qi.hadamard(register[0],1)
    register[1] = qi.hadamard(register[1],1)   

    #Combining the states 
    combined = np.kron(register[0],register[1])


    #Applying the quantum oracle
    post_oracle = np.dot(oracle, combined)


    #calculating the final state
    result = np.dot(np.kron(H,np.eye(2)),post_oracle)

    return result

res = quantumDeutsch()

if res[0]!=0:
    print("Constant")
else:
    print("Balanced")

