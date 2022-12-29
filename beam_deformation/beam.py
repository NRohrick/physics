#Filename: beam.py
#Date: 20221214
#Author: Nathaniel Rohrick
#Description: Models and plots a 2D beam deflection


#Imports
import numpy as np
from scipy.linalg import eigh
#Constants



#Material
class Material():
	def __init__(self, name, temperature, E, v)
		self.name = name
		self.temperature = temperature
		self.E = E
		self.v = v
		self.rho = rho




temperature = 25.0 #degC
E_steel = 200e9
v_steel = 0.25
steel = Material('Steel', temperature, E_steel, v_steel)


class Beam():
	def init__(self, n_elements):
		self.n_elements = n_elements
		#Element Stiffness Matrix
		k = np.array([1,-1],[-1,1]) * (n_elements)
		#Element Mass Matrix
		m = np.array([2,1],[1,2]) / (6*n_elements)

		#Stiffness matrix
		K = np.zeros((n_elements+1, n_elements+1))

		#Mass matrix
		M = np.zeros((n_elements+1, n_elements+1))
		for i in range(n_elements):
			K_temp = np.zeros((n_elements+1, n_elements+1))
			M_temp = np.zeros((n_elements+1, n_elements+1))

			M_temp[i:i+2,i:i+2] = m
			K_temp[i:i+2,i:i+2] = k

			M += M_temp
			K += K_temp

		#Remove the foxed degrees of freedom
		for dof in restrained_dofs:
			for i in [0,1]:
				M = np.delete(M, dof, axis=i)
				K = np.delete(K, dof, axis=i)

		#Eigenvalue problem
		evals, evecs = eigh(K,M)
		frequencies = np.sqrt(evals)
		return M, K, frequencies, evecs