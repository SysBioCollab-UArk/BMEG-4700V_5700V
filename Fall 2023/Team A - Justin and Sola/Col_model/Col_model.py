from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers
Monomer('M1', ['b', 'state'], {'state': ['A', 'I']})
Monomer('M2', ['b', 'state'], {'state': ['A', 'I']})
Monomer('M8', ['b', 'state'], {'state': ['A', 'I']})
Monomer('M9', ['b', 'state'], {'state': ['A', 'I']})
Monomer('T1', ['b', 'state'], {'state': ['A', 'I']})
Monomer('T2', ['b', 'state'], {'state': ['A', 'I']})
Monomer('T4', ['b', 'state'], {'state': ['A', 'I']})

Monomer('C1', ['b', 'state'], {'state': ['A', 'D']})

# Initials
Initial(M1(b=None, state='A'), Parameter('M1_init', 2.3500))
Initial(M2(b=None, state='A'), Parameter('M2_init', 363.12))
Initial(M8(b=None, state='A'), Parameter('M8_init', 2.0800))
Initial(M9(b=None, state='A'), Parameter('M9_init', 85.570))
Initial(T1(b=None, state='A'), Parameter('T1_init', 75.160))
Initial(T2(b=None, state='A'), Parameter('T2_init', 81.740))
Initial(T4(b=None, state='A'), Parameter('T4_init', 0.9700))

Initial(C1(b=None, state='A'), Parameter('C1_init', 44.58))

# Paramaters
# collegen degridation
Parameter('M1C1d', 1)
Parameter('M2C1d', 1)
Parameter('M8C1d', 1)
Parameter('M9C1d', 1)

# Inhibition parameters
Parameter('M1i', 1)
Parameter('M2i', 1)
Parameter('M8i', 1)
Parameter('M9i', 1)
Parameter('T1i', 1)
Parameter('T2i', 1)
Parameter('T4i', 1)

# Forward parameters
Parameter('M1C1f', 1)
Parameter('M2C1f', 1)
Parameter('M8C1f', 1)
Parameter('M9C1f', 1)
Parameter('M1T1f', 1)
Parameter('M2T1f', 1)
Parameter('M8T1f', 1)
Parameter('M9T1f', 1)
Parameter('M1T2f', 1)
Parameter('M2T2f', 1)
Parameter('M8T2f', 1)
Parameter('M9T2f', 1)
Parameter('M1T4f', 1)
Parameter('M2T4f', 1)
Parameter('M8T4f', 1)
Parameter('M9T4f', 1)
# Reverse parameters
Parameter('M1C1r', 1)
Parameter('M2C1r', 1)
Parameter('M8C1r', 1)
Parameter('M9C1r', 1)
Parameter('M1T1r', 1)
Parameter('M2T1r', 1)
Parameter('M8T1r', 1)
Parameter('M9T1r', 1)
Parameter('M1T2r', 1)
Parameter('M2T2r', 1)
Parameter('M8T2r', 1)
Parameter('M9T2r', 1)
Parameter('M1T4r', 1)
Parameter('M2T4r', 1)
Parameter('M8T4r', 1)
Parameter('M9T4r', 1)

# Rules
Rule('M1_binds_C1', M1(b=None, state='A') + C1(b=None, state='A') | M1(b=1, state='A') % C1(b=1, state='A'), M1C1f, M1C1r)
Rule('M2_binds_C1', M2(b=None, state='A') + C1(b=None, state='A') | M2(b=1, state='A') % C1(b=1, state='A'), M2C1f, M2C1r)
Rule('M8_binds_C1', M8(b=None, state='A') + C1(b=None, state='A') | M8(b=1, state='A') % C1(b=1, state='A'), M8C1f, M8C1r)
Rule('M9_binds_C1', M9(b=None, state='A') + C1(b=None, state='A') | M9(b=1, state='A') % C1(b=1, state='A'), M9C1f, M9C1r)

Rule('M1_deg_C1', M1(b=1) % C1(b=1, state='A') >> M1(b=None, state='A') + C1(b=None, state='D'), M1C1d)
Rule('M2_deg_C1', M2(b=1) % C1(b=1, state='A') >> M2(b=None, state='A') + C1(b=None, state='D'), M2C1d)
Rule('M8_deg_C1', M8(b=1) % C1(b=1, state='A') >> M8(b=None, state='A') + C1(b=None, state='D'), M8C1d)
Rule('M9_deg_C1', M9(b=1) % C1(b=1, state='A') >> M9(b=None, state='A') + C1(b=None, state='D'), M9C1d)

Rule('M1_binds_T1', M1(b=None, state='A') + T1(b=None, state='A') | M1(b=1, state='A') % T1(b=1, state='A'), M1T1f, M1T1r)
Rule('M2_binds_T1', M2(b=None, state='A') + T1(b=None, state='A') | M2(b=1, state='A') % T1(b=1, state='A'), M2T1f, M2T1r)
Rule('M8_binds_T1', M8(b=None, state='A') + T1(b=None, state='A') | M8(b=1, state='A') % T1(b=1, state='A'), M8T1f, M8T1r)
Rule('M9_binds_T1', M9(b=None, state='A') + T1(b=None, state='A') | M9(b=1, state='A') % T1(b=1, state='A'), M9T1f, M9T1r)
Rule('M1_binds_T2', M1(b=None, state='A') + T2(b=None, state='A') | M1(b=1, state='A') % T2(b=1, state='A'), M1T2f, M1T2r)
Rule('M2_binds_T2', M2(b=None, state='A') + T2(b=None, state='A') | M2(b=1, state='A') % T2(b=1, state='A'), M2T2f, M2T2r)
Rule('M8_binds_T2', M8(b=None, state='A') + T2(b=None, state='A') | M8(b=1, state='A') % T2(b=1, state='A'), M8T2f, M8T2r)
Rule('M9_binds_T2', M9(b=None, state='A') + T2(b=None, state='A') | M9(b=1, state='A') % T2(b=1, state='A'), M9T2f, M9T2r)
Rule('M1_binds_T4', M1(b=None, state='A') + T4(b=None, state='A') | M1(b=1, state='A') % T4(b=1, state='A'), M1T4f, M1T4r)
Rule('M2_binds_T4', M2(b=None, state='A') + T4(b=None, state='A') | M2(b=1, state='A') % T4(b=1, state='A'), M2T4f, M2T4r)
Rule('M8_binds_T4', M8(b=None, state='A') + T4(b=None, state='A') | M8(b=1, state='A') % T4(b=1, state='A'), M8T4f, M8T4r)
Rule('M9_binds_T4', M9(b=None, state='A') + T4(b=None, state='A') | M9(b=1, state='A') % T4(b=1, state='A'), M9T4f, M9T4r)

Rule('M1_inhib', M1(b=None, state='A') >> M1(b=None, state='I'), M1i)
Rule('M2_inhib', M2(b=None, state='A') >> M2(b=None, state='I'), M2i)
Rule('M8_inhib', M8(b=None, state='A') >> M8(b=None, state='I'), M8i)
Rule('M9_inhib', M9(b=None, state='A') >> M9(b=None, state='I'), M9i)
Rule('T1_inhib', T1(b=None, state='A') >> T1(b=None, state='I'), T1i)
Rule('T2_inhib', T2(b=None, state='A') >> T2(b=None, state='I'), T2i)
Rule('T4_inhib', T4(b=None, state='A') >> T4(b=None, state='I'), T4i)

# Observables
Observable('actC1', C1(b=None, state='A'))
Observable('decC1', C1(b=None, state='D'))

Observable('unbound_M1', M1(b=None, state='A'))
Observable('unbound_M2', M2(b=None, state='A'))
Observable('unbound_M8', M8(b=None, state='A'))
Observable('unbound_M9', M9(b=None, state='A'))
Observable('unbound_T1', T1(b=None, state='A'))
Observable('unbound_T2', T2(b=None, state='A'))
Observable('unbound_T4', T4(b=None, state='A'))

# Simulation commands
tspan = np.linspace(0, 3, 501)
sim = ScipyOdeSimulator(model, tspan, verbose=True)


print(len(model.species))
for i, sp in enumerate (model.species):
    print(i,sp)
print()
for i, ode in enumerate(model.odes):
    print(model.species[i],ode)








output = sim.run()

for obs in model.observables:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc=0)

plt.show()
