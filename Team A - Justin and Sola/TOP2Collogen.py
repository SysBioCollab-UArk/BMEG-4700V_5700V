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
Initial(M1(b=None), Parameter('M1_init', 25))
Initial(M2(b=None), Parameter('M2_init', 25))
Initial(M8(b=None), Parameter('M8_init', 25))
Initial(M9(b=None), Parameter('M9_init', 25))
Initial(T1(b=None), Parameter('T1_init', 25))
Initial(T2(b=None), Parameter('T2_init', 25))
Initial(T4(b=None), Parameter('T4_init', 25))

Initial(C1(b=None, state='A'), Parameter('C1_init', 100))

# Paramaters
# collegen degridation
Parameter('Reactaa', .0001)
Parameter('Reactab', .0001)
Parameter('Reactac', .0001)
Parameter('Reactad', .0001)

# Forward parameters
Parameter('Reactae', .01)
Parameter('Reactag', .01)
Parameter('Reactai', .01)
Parameter('Reactak', .01)
Parameter('Reactam', .01)
Parameter('Reactao', .01)
Parameter('Reactaq', .01)
Parameter('Reactas', .01)
Parameter('Reactau', .01)
Parameter('Reactaw', .01)
Parameter('Reactay', .01)
Parameter('Reactba', .01)
Parameter('Reactbc', .01)
Parameter('Reactbd', .01)
Parameter('Reactbe', .01)
Parameter('Reactbf', .01)
# Reverse parameters
Parameter('Reactaf', .0001)
Parameter('Reactah', .0001)
Parameter('Reactaj', .0001)
Parameter('Reactal', .0001)
Parameter('Reactan', .0001)
Parameter('Reactap', .0001)
Parameter('Reactar', .0001)
Parameter('Reactat', .0001)
Parameter('Reactav', .0001)
Parameter('Reactax', .0001)
Parameter('Reactaz', .0001)
Parameter('Reactbb', .0001)
Parameter('Reactbg', .0001)
Parameter('Reactbh', .0001)
Parameter('Reactbi', .0001)
Parameter('Reactbj', .0001)

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
Observable('unbound_M1', M1(b=None))
Observable('unbound_T1', T1(b=None))

Observable('bound_M1', M1(b=ANY))
Observable('bound_T1', T1(b=ANY))

Observable('decC1', C1(b=None, state='D'))
Observable('actC1', C1(b=None, state='A'))

# Simulation commands
tspan = np.linspace(0, 2900, 501)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
output = sim.run()

for obs in model.observables:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc=0)

plt.show()
