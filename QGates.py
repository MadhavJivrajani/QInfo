import numpy as np 
import math 

"""Pauli Gates X, Y, Z"""
X = np.array([[0,1],[1,0]])
Y = np.array([[0,-1j],[-1j,0]])
Z = np.array([[1,0],[0,-1]])

"""Hadamard Gate H"""
H = (1/np.sqrt(2))*np.array([[1,1],[1,-1]])

"""Phase Gate S"""
S = np.array([[1,0],[0,-1j]])

"""π/8 Gate T"""
T_pi = np.array([[1,0],[0,math.e**((math.pi/4)*1j)]])

"""T dagger"""
T_dag = np.array([[1,0],[0,math.e**((math.pi/4)*-1j)]])

"""Controlled NOT, CNOT"""
cnot = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])

"""SWAP Gate swap"""
swap = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])

# """Controlled Z cZ"""
# cz = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])

"""Controlled Phase Gate cS"""
cs = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1j]])

"""Double Controlled NOT or Toffoli, ccNOT"""
ccnot = np.array([[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,1,0]])

"""Controlled Swap, Fredkin Gate, cSwap"""
cswap = np.array([[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1]])

"""Controlled Pauli Z"""
cz = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])

"""Square root of NOT"""
sqrt_not = np.array([[(1+1j)/2 , (1-1j)/2],[(1-1j)/2 , (1+1j)/2]])

"""P"""
p = np.array([[1,0],[0,1j]])
