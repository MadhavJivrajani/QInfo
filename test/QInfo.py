"""
A library written in Python 3.6.6 for Quantum Information. 
Contains basic functions to perform operations on qubits.
"""

__author__ = "Madhav Jivrajani"

import numpy as np 
import math
from QGates import *
from mapping import mapping_2 as m2
from mapping import mapping_2_rev as m2_rev

class Qubit:
    def __init__(self, n):
        """Initializes a qubit in the |0> state"""
        self.n = n
        self.state = self.construct_standard_basis(1)[0]
        self.measureRes = {}

    def initQubit(self):
        """Initializes qubit to |0> state"""
        return self.state

    def to_bra(self, ket):
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
                self.state = (ket.conjugate()).transpose()
            except:
                print("Invalid type: Argument passed must be an object of type numpy.ndarray.")


    def to_ket(self, bra):
        """
        Converts a bra to the corresponding ket.
        bra should be aan object of type numpy.ndarray
        returns the ket which is also an object of type numpy.ndarray.
        """
        if isinstance(bra, np.ndarray):
                self.state = (bra.conjugate()).transpose()
        else:
            try:
                bra = np.array(bra)
                self.state = (bra.conjugate()).transpose()
            except:
                print("Invalid type: Argument passed must be an object of type numpy.ndarray.")

    def inner_product(self, psi_1, psi_2):
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

    def check_validity(self, psi):
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

    def construct_standard_basis(self, n):
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

    def normalising_factor(self, n):
        """
        Returns the factor by which a state is to be multiplied to normalise it.
        """
        if n>=0:
            return 1/math.sqrt(2**n)
        else:
            print("Number of Qubits must be a non-negative number.")

    def equal_superposition(self, states):
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
        
        norm_factor = np.array(self.normalising_factor(math.log(len(states))/math.log(2))) #since, len = 2^n
        superposition = states[0]*norm_factor
        for state in states[1:]:
            superposition = superposition + state * norm_factor
        return superposition

    def tensor_combine(self, states):
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

    def measure_standard(self):
        """
        Measures a qubit in the standard basis. 
        Returns a dictionary with probability values. 
        Keys are 0, 1 indicative of the standard basis and values are the corresponding probability values.
        psi must be an object of type numpy.ndarray
        """
        if isinstance(self.state, np.ndarray):
            basis = self.construct_standard_basis(1)
            p0 = float(self.inner_product(self.state, basis[0])[0])
            p1 = float(self.inner_product(self.state, basis[1])[0])
            self.measureRes = {"0" : np.square(np.abs(p0)), "1" : np.square(np.abs(p1))}    
        else:
            try:
                self.state = np.array(self.state)
                basis = self.construct_standard_basis(1)
                p0 = float(self.inner_product(self.state, basis[0])[0])
                p1 = float(self.inner_product(self.state, basis[1])[0])
                self.measureRes = {"0" : np.square(np.abs(p0)), "1" : np.square(np.abs(p1))}    
            except:
                print("Invalid type: Argument passed must be an object of type numpy.ndarray.")

    def measure(self, state):
        """
        Measures the state of psi in another state.
        Both psi and state must be objects of type numpy.ndarray.
        Returns a dictionary with key as "state" and value as the probability.
        """
        if isinstance(self.state, np.ndarray) and isinstance(state, np.ndarray):
            p_state = self.inner_product(self.state, state)[0][0]
            return {"state" : np.square(np.abs(p_state))}
        else:
            try:
                self.state = np.array(self.state)
                state = np.array(state)
                p_state = self.inner_product(self.state, state)[0][0]
                return {"state" : np.square(np.abs(p_state))}
            except:
                print("Invalid type: Arguments passed must be objects of type numpy.ndarray.")

    def hadamard(self):
        """
        Application of 1 qubit hadamard transform.
        """
        if not isinstance(self.state, np.ndarray):
            print("Invalid type: state must be object of type numpy.ndarray.")
            return
        
        self.state = np.dot(H, self.state)

    def pauliX(self):
        """
        Application of 1 qubit Pauli X transform.
        """
        if not isinstance(self.state, np.ndarray):
            print("Invalid type: state must be object of type numpy.ndarray.")
            return

        self.state = np.dot(X, self.state)   

    def pauliY(self):
        """
        Application of 1 qubit Pauli Y transform.
        """
        if not isinstance(self.state, np.ndarray):
            print("Invalid type: state must be object of type numpy.ndarray.")
            return

        self.state = np.dot(Y, self.state)   

    def pauliZ(self):
        """
        Application of 1 qubit Pauli Z transform.
        """
        if not isinstance(self.state, np.ndarray):
            print("Invalid type: state must be object of type numpy.ndarray.")
            return

        self.state = np.dot(Z, self.state)   

    def phase(self):
        """
        Application phase measurement on 1 qubit.
        """
        if not isinstance(self.state, np.ndarray):
            print("Invalid type: state must be object of type numpy.ndarray.")
            return

        self.state = np.dot(S, self.state)   

    def T(self):
        """
        Application 1 qubit π/8 or T gate.
        """
        if not isinstance(self.state, np.ndarray):
            print("Invalid type: state must be object of type numpy.ndarray.")
            return

        self.state = np.dot(T_pi, self.state)   

    def cNOT(self, this):
        """
        Application of a controlled not.
        Takes a Qubit object as argument which serves the purpose of the controlled qubit.
        """
        if not isinstance(self.state, np.ndarray):
            self.state = np.array(self.state)
        if not isinstance(this.state, np.ndarray):
            this.state = np.array(this.state)
        if m2[self.state.tobytes()]==1:
            self.pauliX()

    def cZ(self, this):
        """
        Application of a controlled Pauli Z.
        Takes a Qubit object as argument which serves the purpose of the controlled qubit.
        """
        if not isinstance(self.state, np.ndarray):
            self.state = np.array(self.state)
        if not isinstance(this.state, np.ndarray):
            this.state = np.array(this.state)
        if m2[self.state.tobytes()]==1:
            self.pauliZ()
    
    def cY(self, this):
        """
        Application of a controlled Pauli Y.
        Takes a Qubit object as argument which serves the purpose of the controlled qubit.
        """
        if not isinstance(self.state, np.ndarray):
            self.state = np.array(self.state)
        if not isinstance(this.state, np.ndarray):
            this.state = np.array(this.state)
        if m2[self.state.tobytes()]==1:
            self.pauliY()

    def ccNOT(self, this1, this2):
        """
        Application of a double controlled not.
        this1, this2 Qubit objects act as controls 
        self is target Qubit object.
        """
        if not isinstance(self.state, np.ndarray):
            self.state = np.array(self.state)
        if not isinstance(this1.state, np.ndarray):
            this1.state = np.array(this1.state)
        if not isinstance(this2.state, np.ndarray):
            this2.state = np.array(this2.state)
        
        if m2[self.state.tobytes()]==1 and m2[this1.state.tobytes()]==1:
            self.pauliX()

    def combineGates(self, gates):
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

    def applyGate(self, state, gate):
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

    def densityMatrix(self, psi):
        """
        Forms the density matrix representation of a qubit.
        """
        if not isinstance(psi, np.ndarray):
            psi = np.array(psi)
        else:
            psi_bra = self.to_bra(psi)
            return psi * psi_bra

    def measureDensity(self, psi, state):
        """
        Measurement using density matrices
        """
        if not isinstance(psi, np.ndarray):
            psi = np.array(psi)
        if not isinstance(state, np.ndarray):
            state = np.array(state)
        state_bra = self.to_bra(state)
        rho = self.densityMatrix(psi)
        return {"state" : np.trace(state_bra * rho * state)}

class Register:
    def __init__(self, n):
        """Initializes a register with all qubits in |0> state"""
        self.n = n
        self.register = np.array([Qubit(self.n) for _ in range(self.n)])
        

    def initialize(self):
        """Returns Quantum Register"""
        
        return self.register
