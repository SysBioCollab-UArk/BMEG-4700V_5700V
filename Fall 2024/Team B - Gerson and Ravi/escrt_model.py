from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers
Monomer('cargo', ['escrt0', 'state'], {'state': ['x', 'ub']})
Monomer('ESCRT_0', ['cargo', 'escrt1'])
Monomer('ESCRT_I', ['escrt0', 'escrt2'])
Monomer('ESCRT_II', ['escrt1', 'escrt2'])
# Monomer('ESCRT_III')

# Initial conditions
Parameter('cargo_init', 1e5)
Parameter('ESCRT_0_init', 100)
Parameter('ESCRT_I_init', 100)
Parameter('ESCRT_II_init', 100)

Initial(cargo(escrt0=None, state='x'), cargo_init)
Initial(ESCRT_0(cargo=None, escrt1=None), ESCRT_0_init)
Initial(ESCRT_I(escrt0=None, escrt2=None), ESCRT_I_init)
Initial(ESCRT_II(escrt1=None, escrt2=None), ESCRT_II_init)

# Rules
Parameter('k0', 1)
Rule('cargo_ubiquitinated',
     cargo(escrt0=None, state='x') >> cargo(escrt0=None, state='ub'), k0)

Parameter('kf1', 1)
Parameter('kr1', 1)
Rule('ESCRT0_binds_ubiq_cargo',
     ESCRT_0(cargo=None, escrt1=None) + cargo(escrt0=None, state='ub') |
     ESCRT_0(cargo=1, escrt1=None) % cargo(escrt0=1, state='ub'), kf1, kr1)

Parameter('kf2', 1)
Parameter('kr2', 1)
Rule('ESCRT0_recruits_ESCRT1',
     ESCRT_0(cargo=1, escrt1=None) % cargo(escrt0=1, state='ub')
     + ESCRT_I(escrt0=None, escrt2=None) |
     ESCRT_0(cargo=1, escrt1=2) % cargo(escrt0=1, state='ub')
     % ESCRT_I(escrt0=2, escrt2=None), kf2, kr2)

Parameter('kf3', 1)
Parameter('kr3', 1)
Rule('ESCRT1_recruits_ESCRT2',
     ESCRT_0(cargo=1, escrt1=2) % cargo(escrt0=1, state='ub')
     % ESCRT_I(escrt0=2, escrt2=None) + ESCRT_II(escrt1=None, escrt2=None) |
     ESCRT_0(cargo=1, escrt1=2) % cargo(escrt0=1, state='ub')
     % ESCRT_I(escrt0=2, escrt2=3) % ESCRT_II(escrt1=3, escrt2=None), kf3, kr3)

Parameter('kf4', 1)
Parameter('kr4', 1)
Rule('membrane_invagination',
     ESCRT_II(escrt1=ANY, escrt2=None) + ESCRT_II(escrt1=ANY, escrt2=None) |
     ESCRT_II(escrt1=ANY, escrt2=1) % ESCRT_II(escrt1=ANY, escrt2=1), kf4, kr4)

# Rule('..', ...)
# ...
# ...

# Observables
# Observable('cargo_tot', cargo())
# Observable('cargo_ub', cargo(state='ub'))
Observable('cargo_ESCRT0', cargo(escrt0=ANY))
Observable('cargo_ESCRT0_ESCRT1', cargo(escrt0=1) % ESCRT_0(cargo=1, escrt1=ANY))
Observable('cargo_ESCRT0_ESCRT1_ESCRT2',
           cargo(escrt0=1) % ESCRT_0(cargo=1, escrt1=2) % ESCRT_I(escrt0=2, escrt2=ANY))
Observable('initial_bud', ESCRT_II(escrt2=ANY))

# simulation commands
tspan = np.linspace(0, 1, 101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
output = sim.run()

for obs in model.observables:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc='best')
plt.tight_layout()

plt.show()












