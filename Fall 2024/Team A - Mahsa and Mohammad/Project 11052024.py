from pysb import *
from pysb.examples.explicit import model
from pysb.simulator import scipyode
import numpy as np
import matplotlib.pyplot as plt
from sympy.stats.rv import sample_iter_subs

Model()

# Monomers
Monomer(name='A', sites=['b'])
Monomer(name='B', sites=['a','state'], site_states:{'states':['u','p']})




# Initial conditions


# Rules



# Observables


# Simulation commands

