from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers

Monomer('A', ['b'])
Monomer('B', ['a','state'], {'state':['u','p']})

# Initial Conditions

Parameter('A_0', 100)
Parameter('B_0', 80)
Initial(A(b=None), A_0)
Initial(B(a=None, state='u'), B_0)

# Rules

Parameter('kf_AB', 1)
Parameter('kr_AB', 10)
Parameter('k_B_phos', 0.1)
Parameter(' k_ABp_unbind', 1)
Parameter( 'k_Bp_unphos', 10)
Rule('A_B_binds', A(b=None), B(a=None,state='u') | A(b=1) % B(a=1, state='u'), kf_AB, kr_AB)
Rule('AB_phos', A(b=1) % B(a=1, state='u') >> A(b=1) % B(a=1, state='p'), k_B_phos)
Rule('AB_unbinds', A(b=1) % B(a=1, state='p') >> A(b=None) + B(a=None, state='p'), k_ABp_unbind)
Rule('Bp_unphos', B(a=None, state = 'p') >> B(a=None, state = 'u'), k_Bp_unphos)

# Observables




# Simulation Commands
