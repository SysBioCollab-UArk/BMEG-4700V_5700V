from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers
Monomer('E', ['b'])
Monomer('X', ['b', 'state'], {'state': ['S', 'P']})

# Initial Conditions
Initial(E(b=None), Parameter('E_0', 10))
Initial(X(b=None, state='S'), Parameter('S_0', 100))

# Rules
Parameter('kf', 0.1)  # 0.1
Parameter('kr', 1000)
Rule('E_binds_S', E(b=None) + X(b=None, state='S') | E(b=1) % X(b=1, state='S'), kf, kr)

Parameter('kcat', 100)
Rule('ES_makes_P', E(b=1) % X(b=1, state='S') >> E(b=None) + X(b=None, state='P'), kcat)

# Observables
Observable('E_tot', E())
Observable('E_free', E(b=None))
Observable('S_free', X(b=None, state='S'))
Observable('ES_complex', E(b=ANY))
Observable('Product', X(state='P'))
Observable('Better_be_zero', X(b=ANY, state='P'))

# Simulation commands
tspan = np.linspace(0, 50, 501)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
output = sim.run()

for obs in model.observables:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc=0)

plt.show()
