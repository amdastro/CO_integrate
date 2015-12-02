import physical as ph
import numpy as np

'''Formation rate from Lazzati 2015'''
def formation(T):
	K = 4.467e-17/((T/4467.)**(-2.08) + (T/4467.)**(-0.22))**0.5
	return K

'''Destructiom rate by thermal processes'''
def thermal(T):
	K = 4.467e-17/((T/4467.)**(-2.08) + (T/4467.)**(-0.22))**0.5\
		*(ph.h**2 / (2.*np.pi*ph.mu*ph.kB*T))**(-1.5) * np.exp(-ph.B_CO/(ph.kB*T))
	return K

