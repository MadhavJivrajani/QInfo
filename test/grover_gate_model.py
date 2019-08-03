"""
Grover's algorithm using gates, implemented for 2 qubits.
Marked element is 3. 
"""

__author__ = "Madhav Jivrajani"

import QInfo as qi 
import QGates as gates 
import numpy as np 
import math 

reg = [qi.construct_standard_basis(1)[0] for _ in range(2)]
I = np.eye(2)

reg.append(qi.construct_standard_basis(1)[1])
reg[0] = qi.hadamard(reg[0],1)
reg[1] = qi.hadamard(reg[1],1)
reg[2] = qi.hadamard(reg[2],1)

#reg[0] = qi.pauliX(reg[0],1)
#reg[1] = qi.pauliX(reg[1],1)

psi = qi.ccNOT(qi.tensor_combine(np.array(reg)))

gate = qi.combineGates([gates.H,gates.H,np.eye(2)])
psi = qi.applyGate(psi, gate)

gate = qi.combineGates([gates.X,gates.X,np.eye(2)])
psi = qi.applyGate(psi, gate)

gate = qi.combineGates([gates.cz,np.eye(2)])
psi = qi.applyGate(psi, gate)

gate = qi.combineGates([gates.X,gates.X,np.eye(2)])
psi = qi.applyGate(psi, gate)

gate = qi.combineGates([gates.H,gates.H,np.eye(2)])
psi = qi.applyGate(psi, gate)


print(psi)
print(qi.measure(psi, qi.construct_standard_basis(3)[7]))

#TODO : get amplitudes from final state psi.