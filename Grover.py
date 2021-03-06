"""
An implementation of Grover's Algorithm written using the QInfo Python Library. 
"""

__author__ = "Madhav Jivrajani"

from QInfo import Register
from QGates import H
import numpy as np 
import math 
import matplotlib.pyplot as plt 

class Grover:
    def __init__(self, n, target):
        """
        n      :     signifies the number of qubits being used, must be an integer. 
        target :     signifies the element being searched for, for which the amplitude must be amplified. 
        reg    :     signifies a quantum register holding n qubits in the |0> state.
        vec    :     signifies the combined state of all qubits. 
        states :     represents a list storing the state of the combined qubit state at various stages
                     of the algorithm.
        """
        self.n = n
        self.target = target
        self.reg = Register(self.n).initialize()
        self.vec = self.reg[0].tensor_combine(self.reg[1:])
        self.states = []
    
    def oracle(self):
        """Implements the quantum oracle"""
        matZero = np.zeros((2**self.n,2**self.n))
        for i in range(0,2**self.n):
            if i == self.target:
                f = 1
            else:
                f = 0

            matZero[i,i] = (-1)**f
        return matZero
    
    def invertAboutAverage(self):
        """Operator to invert amplitudes about the average."""
        A = np.ones(2**self.n)/(2**self.n)
        invertMatrix = 2*A - np.identity(2**self.n)
        return invertMatrix

    def groverIteration(self):
        """Returns the maximum grover's iterations. Quadratic speedup."""
        return (np.pi/4)*math.sqrt(2**self.n)
    
    def grover(self):
        """Implements the grover's diffusion operator and the grover's algorithm."""
        H_temp = 1
        if self.n == 1:
            H_temp = H
        for _ in range(self.n):
            H_temp = np.kron(H_temp, H)
        self.vec = np.dot(H_temp, self.vec)
        self.states.append(self.vec)
        self.iter = int(np.trunc(self.groverIteration()))
        for _ in range(1,self.iter):
            oracle = self.oracle()
            self.vec = np.dot(oracle, self.vec)
            self.states.append(self.vec)
            invertMatrix = self.invertAboutAverage()
            self.vec = np.dot(invertMatrix, self.vec)
            self.states.append(self.vec)
    
    def plotAlgo(self):
        """Plots bar graphs for probablity amplitudes at various stages of the algorithm."""
        plot_total = len(self.states)
        _, ax = plt.subplots(nrows = plot_total, ncols = 1)
        i = 0
        for row in ax:
            curr_plot = np.array([k[0] for k in self.states[i]])
            index = np.arange(len(curr_plot))
            row.bar(index, curr_plot)
            row.set_xticks(list(range(2**self.n)))
            row.set_xlabel("States")
            row.set_ylabel("Amplitudes")
            i+=1
        plt.show()

    def printAlgo(self):
        """Prints the states and oracles at different stages of the algorithm."""
        print(self.states[0])
        print()
        print("Oracle: \n")
        print(self.oracle())
        for state in self.states[1:]:
            print()
            print(state)
            print()

