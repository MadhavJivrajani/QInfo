"""
A library written in Python 3.6.6 for Quantum Information. 
Contains basic functions to perform operations on qubits.
"""

__author__ = "Madhav Jivrajani"

import numpy as np 
import math
from QGates import *


def to_bra(ket):
    """
    Converts a ket to the corresponding bra.
    ket should be an object of type numpy.ndarray
    returns the bra which is also a numpy array.
    """
    if isinstance(ket, np.ndarray):
        return (ket.conjugate()).transpose()
    else:
        try:
            ket = np.array(ket)
            return (ket.conjugate()).transpose()
        except:
            print("Invalid type: Argument passed must be an object of type numpy.ndarray.")


def to_ket(bra):
    """
    Converts a bra to the corresponding ket.
    bra should be aan object of type numpy.ndarray
    returns the ket which is also an object of type numpy.ndarray.
    """
    if isinstance(bra, np.ndarray):
            return (bra.conjugate()).transpose()
    else:
        try:
            bra = np.array(bra)
            return (bra.conjugate()).transpose()
        except:
            print("Invalid type: Argument passed must be an object of type numpy.ndarray.")

def inner_product(psi_1, psi_2):
    """
    Calculates the inner product of psi_1 with psi_2.
    Inner product can be comprehended as the dot product between two vectors.
    psi_1 and psi_2 must be objects of type numpy.ndarray. 
    Returns a numeric value.
    """
    if isinstance(psi_1, np.ndarray) and isinstance(psi_2, np.ndarray):
        return np.dot((psi_1.conjugate()).transpose(), psi_2)
    else:
        try:
            psi_1 = np.array(psi_1)
            psi_2 = np.array(psi_2)
            return np.dot((psi_1.conjugate()).transpose(), psi_2)
        except:
            print("Invalid type: Arguments passed must be objects of type numpy.ndarray.")

def check_validity(psi):
    """
    Checks whether a vector can be a qubit or not.
    Argument passed must be an object of type numpy.ndarray
    Returns a boolean.
    """
    if isinstance(psi, np.ndarray):
            return np.linalg.norm(psi, 2) == 1
    else:
        try:
            psi = np.array(psi)
            return np.linalg.norm(psi, 2) == 1
        except:
            print("Invalid type: Argument passed must be an object of type numpy.ndarray.")

def construct_standard_basis(n):
    """
    Constructs the standard basis vectors for n qubits. 
    Returns a numpy array with containing the standard basis vectors, each being a numpy array of shape (1,2^n)
    """
    basis = []
    for i in range(2**n):
        a = []
        for j in range(2**n):
            if(j==i):
                a.append([1])
            else:
                a.append([0])
        basis.append(a)
    return np.array(basis)

def normalising_factor(n):
    """
    Returns the factor by which a state is to be multiplied to normalise it.
    """
    if n>=0:
        return 1/math.sqrt(2**n)
    else:
        print("Number of Qubits must be a non-negative number.")

def equal_superposition(states):
    """
    Creates an equal superposition from an array of states.
    states must be an object of type numpy.ndarray.
    each element of states must also be an object of type numpy.ndarray.
    """
    flag = True
    if not isinstance(states, np.ndarray):
        try:
            states = np.array(states)
        except:
            print("Invalid type: Argument passed must be an object of type numpy.ndarray.")
    
    i = 0
    while i<len(states) and flag:
        if not isinstance(states[i], np.ndarray):
            flag = False
        i+=1
    if flag == False:
        try:
            while i<len(states):
                states[i] = np.array(states[i])
        except:
            print("Invalid type: Elements must be objects of type numpy.ndarray.")
    
    norm_factor = np.array(normalising_factor(math.log(len(states))/math.log(2))) #since, len = 2^n
    superposition = states[0]*norm_factor
    for state in states[1:]:
        superposition = superposition + state * norm_factor
    return superposition

def tensor_combine(states):
    """
    Combines multiple qubits using the tensor product. 
    states must be an object of type numpy.ndarray.
    Each element of states must also be an object of type numpy.ndarray.
    Returns the combined state of qubits passed as an array to the function, which is also an object of type numpy.ndarray
    """
    if not len(states)>=2:
        print("Need atleast two states to combine.")
        return
    combined = np.kron(states[0],states[1])
    if len(states)>2:
        for state in states[2:]:
            combined = np.kron(combined, state)
        return combined
    return combined

def measure_standard(psi):
    """
    Measures a qubit in the standard basis. 
    Returns a dictionary with probability values. 
    Keys are 0, 1 indicative of the standard basis and values are the corresponding probability values.
    psi must be an object of type numpy.ndarray
    """
    if isinstance(psi, np.ndarray):
        basis = construct_standard_basis(1)
        p0 = float(inner_product(psi, basis[0])[0])
        p1 = float(inner_product(psi, basis[1])[0])
        return {"0" : np.square(np.abs(p0)), "1" : np.square(np.abs(p1))}    
    else:
        try:
            psi = np.array(psi)
            basis = construct_standard_basis(1)
            p0 = float(inner_product(psi, basis[0])[0])
            p1 = float(inner_product(psi, basis[1])[0])
            return {"0" : np.square(np.abs(p0)), "1" : np.square(np.abs(p1))}    
        except:
            print("Invalid type: Argument passed must be an object of type numpy.ndarray.")

def measure(psi, state):
    """
    Measures the state of psi in another state.
    Both psi and state must be objects of type numpy.ndarray.
    Returns a dictionary with key as "state" and value as the probability.
    """
    if isinstance(psi, np.ndarray) and isinstance(state, np.ndarray):
        p_state = inner_product(psi, state)[0][0]
        return {"state" : np.square(np.abs(p_state))}
    else:
        try:
            psi = np.array(psi)
            state = np.array(state)
            p_state = inner_product(psi, state)[0][0]
            return {"state" : np.square(np.abs(p_state))}
        except:
            print("Invalid type: Arguments passed must be objects of type numpy.ndarray.")

def hadamard(state, n):
    """
    Application of n qubit hadamard transform on an n qubit state.
    state must be an object of type numpy.ndarray
    n must be an integer.
    """
    if not isinstance(state, np.ndarray):
        print("Invalid type: state must be object of type numpy.ndarray.")
        return
    if not type(n)==int:
        print("Invalid type: number of qubits must be an integer.")
        return
    if not math.log(state.shape[0])/math.log(2)==n:
        print("Inconsistent dimension: number of qubits in state does not match number of qubits passed.")
        return        

    if n==1:
        H_temp = H
    else:
        H_temp = 1
        for i in range(0,n):
            H_temp = np.kron(H_temp, H) #kron product of matrices 
    return np.dot(H_temp, state)

def pauliX(state, n):
    """
    Application of n qubit Pauli X transform on an n qubit state.
    state must be an object of type numpy.ndarray
    n must be an integer.
    """
    if not isinstance(state, np.ndarray):
        print("Invalid type: state must be object of type numpy.ndarray.")
        return
    if not type(n)==int:
        print("Invalid type: number of qubits must be an integer.")
        return
    if not math.log(state.shape[0])/math.log(2)==n:
        print("Inconsistent dimension: number of qubits in state does not match number of qubits passed.")
        return        

    if n==1:
        X_temp = X
    else:
        X_temp = 1
        for i in range(0,n):
            X_temp = np.kron(X_temp, X) #kron product of matrices 
    return np.dot(X_temp, state)   

def pauliY(state, n):
    """
    Application of n qubit Pauli Y transform on an n qubit state.
    state must be an object of type numpy.ndarray
    n must be an integer.
    """
    if not isinstance(state, np.ndarray):
        print("Invalid type: state must be object of type numpy.ndarray.")
        return
    if not type(n)==int:
        print("Invalid type: number of qubits must be an integer.")
        return
    if not math.log(state.shape[0])/math.log(2)==n:
        print("Inconsistent dimension: number of qubits in state does not match number of qubits passed.")
        return        

    if n==1:
        Y_temp = Y
    else:
        Y_temp = 1
        for i in range(0,n):
            Y_temp = np.kron(Y_temp, Y) #kron product of matrices 
    return np.dot(Y_temp, state)

def pauliZ(state, n):
    """
    Application of n qubit Pauli Z transform on an n qubit state.
    state must be an object of type numpy.ndarray
    n must be an integer.
    """
    if not isinstance(state, np.ndarray):
        print("Invalid type: state must be object of type numpy.ndarray.")
        return
    if not type(n)==int:
        print("Invalid type: number of qubits must be an integer.")
        return
    if not math.log(state.shape[0])/math.log(2)==n:
        print("Inconsistent dimension: number of qubits in state does not match number of qubits passed.")
        return        

    if n==1:
        Z_temp = Z
    else:
        Z_temp = 1
        for i in range(0,n):
            Z_temp = np.kron(Z_temp, Z) #kron product of matrices 
    return np.dot(Z_temp, state)

def phase(state, n):
    """
    Application n qubit phase measurement.
    state must be an object of type numpy.ndarray
    n must be an integer.
    """
    if not isinstance(state, np.ndarray):
        print("Invalid type: state must be object of type numpy.ndarray.")
        return
    if not type(n)==int:
        print("Invalid type: number of qubits must be an integer.")
        return
    if not math.log(state.shape[0])/math.log(2)==n:
        print("Inconsistent dimension: number of qubits in state does not match number of qubits passed.")
        return        

    if n==1:
        S_temp = S
    else:
        S_temp = 1
        for i in range(0,n):
            S_temp = np.kron(S_temp, S) #kron product of matrices 
    return np.dot(S_temp, state)

def T(state, n):
    """
    Application n qubit π/8 or T gate.
    state must be an object of type numpy.ndarray
    n must be an integer.
    """
    if not isinstance(state, np.ndarray):
        print("Invalid type: state must be object of type numpy.ndarray.")
        return
    if not type(n)==int:
        print("Invalid type: number of qubits must be an integer.")
        return
    if not math.log(state.shape[0])/math.log(2)==n:
        print("Inconsistent dimension: number of qubits in state does not match number of qubits passed.")
        return        

    if n==1:
        T_temp = T_pi
    else:
        T_temp = 1
        for i in range(0,n):
            T_temp = np.kron(T_temp, T_pi) #kron product of matrices 
    return np.dot(T_temp, state)

def cNOT(combined):
    """
    Application of a controlled not.
    combined must be object of type numpy.ndarray.
    combined is the combined state of the controlling and controlled qubit states. 
    combined = controlling ⊗ controlled.
    """
    if isinstance(combined, np.ndarray):
        try:
            return np.dot(cnot, combined)
        except Exception as e:
            raise e

    else:
        try:
            combined = np.array(combined)
            try:
                return np.dot(cnot, combined)
            except Exception as e:
                raise e

        except:
            print("Invalid type: Argument passed must be an object of type numpy.ndarray.")

def cZ(combined):
    """
    Application of a controlled Z.
    combined must be object of type numpy.ndarray.
    combined is the combined state of the controlling and controlled qubit states. 
    combined = controlling ⊗ controlled.
    """
    if isinstance(combined, np.ndarray):
        try:
            return np.dot(cz, combined)
        except Exception as e:
            raise e

    else:
        try:
            combined = np.array(combined)
            try:
                return np.dot(cz, combined)
            except Exception as e:
                raise e

        except:
            print("Invalid type: Argument passed must be an object of type numpy.ndarray.")

def ccNOT(combined):
    """
    Application of a double controlled not.
    combined must be object of type numpy.ndarray.
    combined is the combined state of the controlling and controlled qubit states. 
    combined = controlling ⊗ controlled.
    """
    if isinstance(combined, np.ndarray):
        try:
            return np.dot(ccnot, combined)
        except Exception as e:
            raise e

    else:
        try:
            combined = np.array(combined)
            try:
                return np.dot(ccnot, combined)
            except Exception as e:
                raise e

        except:
            print("Invalid type: Argument passed must be an object of type numpy.ndarray.")

def combineGates(gates):
    """
    Combines gates in the same step of a quantum circuit. 
    If not gate is present, identity matrix of dimension n x n to be passed,
    where n is the number of qubits in that wire.
    """
    if not len(gates)>=2:
        print("Need atleast two gates to combine.")
        return
    combined = np.kron(gates[0],gates[1])
    if len(gates)>2:
        for gate in gates[2:]:
            combined = np.kron(combined, gate)
        return combined
    return combined

def applyGate(state, gate):
    """
    Applies gate to state.
    Both gate and state must be objects of type numpy.ndarray
    Returns an object of type numpy.ndarray representing the combined state. 
    """
    if not isinstance(state, np.ndarray):
        state = np.array(state)
    if not isinstance(gate, np.ndarray):
        gate = np.array(gate)
    return np.dot(gate, state)