from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers
Monomer('M1', ['b'])
Monomer('M2', ['b'])
Monomer('M8', ['b'])
Monomer('M9', ['b'])
Monomer('T1', ['b'])
Monomer('T2', ['b'])
Monomer('T4', ['b'])

Monomer('C1', ['b', 'state'], {'state': ['A', 'D']})

# Initials
Initial(M1(b=None), Parameter('M1_init', 2.3500))
Initial(M2(b=None), Parameter('M2_init', 363.12))
Initial(M8(b=None), Parameter('M8_init', 2.0800))
Initial(M9(b=None), Parameter('M9_init', 85.570))
Initial(T1(b=None), Parameter('T1_init', 75.160))
Initial(T2(b=None), Parameter('T2_init', 81.740))
Initial(T4(b=None), Parameter('T4_init', 0.9700))

Initial(C1(b=None, state='A'), Parameter('C1_init', 44.58))

# Paramaters
# collegen degridation
Parameter('Reactaa', 1)
Parameter('Reactab', 1)
Parameter('Reactac', 1)
Parameter('Reactad', 1)

# Forward parameters
Parameter('Reactae', 1)
Parameter('Reactag', 1)
Parameter('Reactai', 1)
Parameter('Reactak', 1)
Parameter('Reactam', 1)
Parameter('Reactao', 1)
Parameter('Reactaq', 1)
Parameter('Reactas', 1)
Parameter('Reactau', 1)
Parameter('Reactaw', 1)
Parameter('Reactay', 1)
Parameter('Reactba', 1)
Parameter('Reactbc', 1)
Parameter('Reactbd', 1)
Parameter('Reactbe', 1)
Parameter('Reactbf', 1)
# Reverse parameters
Parameter('Reactaf', 1)
Parameter('Reactah', 1)
Parameter('Reactaj', 1)
Parameter('Reactal', 1)
Parameter('Reactan', 1)
Parameter('Reactap', 1)
Parameter('Reactar', 1)
Parameter('Reactat', 1)
Parameter('Reactav', 1)
Parameter('Reactax', 1)
Parameter('Reactaz', 1)
Parameter('Reactbb', 1)
Parameter('Reactbg', 1)
Parameter('Reactbh', 1)
Parameter('Reactbi', 1)
Parameter('Reactbj', 1)

# Rules
Rule('M1_binds_C1', M1(b=None) + C1(b=None, state='A') | M1(b=1) % C1(b=1, state='A'), Reactbc, Reactbg)
Rule('M2_binds_C1', M2(b=None) + C1(b=None, state='A') | M2(b=1) % C1(b=1, state='A'), Reactbd, Reactbh)
Rule('M8_binds_C1', M8(b=None) + C1(b=None, state='A') | M8(b=1) % C1(b=1, state='A'), Reactbe, Reactbi)
Rule('M9_binds_C1', M9(b=None) + C1(b=None, state='A') | M9(b=1) % C1(b=1, state='A'), Reactbf, Reactbj)

Rule('M1_deg_C1', M1(b=1) % C1(b=1, state='A') >> M1(b=None) + C1(b=None, state='D'), Reactaa)
Rule('M2_deg_C1', M2(b=1) % C1(b=1, state='A') >> M2(b=None) + C1(b=None, state='D'), Reactab)
Rule('M8_deg_C1', M8(b=1) % C1(b=1, state='A') >> M8(b=None) + C1(b=None, state='D'), Reactac)
Rule('M9_deg_C1', M9(b=1) % C1(b=1, state='A') >> M9(b=None) + C1(b=None, state='D'), Reactad)

Rule('M1_binds_T1', M1(b=None) + T1(b=None) | M1(b=1) % T1(b=1), Reactae, Reactaf)
Rule('M2_binds_T1', M2(b=None) + T1(b=None) | M2(b=1) % T1(b=1), Reactag, Reactah)
Rule('M8_binds_T1', M8(b=None) + T1(b=None) | M8(b=1) % T1(b=1), Reactai, Reactaj)
Rule('M9_binds_T1', M9(b=None) + T1(b=None) | M9(b=1) % T1(b=1), Reactak, Reactal)
Rule('M1_binds_T2', M1(b=None) + T2(b=None) | M1(b=1) % T2(b=1), Reactam, Reactan)
Rule('M2_binds_T2', M2(b=None) + T2(b=None) | M2(b=1) % T2(b=1), Reactao, Reactap)
Rule('M8_binds_T2', M8(b=None) + T2(b=None) | M8(b=1) % T2(b=1), Reactaq, Reactar)
Rule('M9_binds_T2', M9(b=None) + T2(b=None) | M9(b=1) % T2(b=1), Reactas, Reactat)
Rule('M1_binds_T4', M1(b=None) + T4(b=None) | M1(b=1) % T4(b=1), Reactau, Reactav)
Rule('M2_binds_T4', M2(b=None) + T4(b=None) | M2(b=1) % T4(b=1), Reactaw, Reactax)
Rule('M8_binds_T4', M8(b=None) + T4(b=None) | M8(b=1) % T4(b=1), Reactay, Reactaz)
Rule('M9_binds_T4', M9(b=None) + T4(b=None) | M9(b=1) % T4(b=1), Reactba, Reactbb)

# Observables
Observable('actC1', C1(b=None, state='A'))
Observable('decC1', C1(b=None, state='D'))

Observable('unbound_M1', M1(b=None))
Observable('unbound_M2', M2(b=None))
Observable('unbound_M8', M8(b=None))
Observable('unbound_M9', M9(b=None))
Observable('unbound_T1', T1(b=None))
Observable('unbound_T2', T2(b=None))
Observable('unbound_T4', T4(b=None))
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
