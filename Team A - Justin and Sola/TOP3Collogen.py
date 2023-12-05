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
Parameter('Reactaa', 1)
Parameter('Reactab', 1)
Parameter('Reactac', 1)
Parameter('Reactad', 1)

# Inhibition parameters
Parameter('Reactbk', 1)
Parameter('Reactbl', 1)
Parameter('Reactbm', 1)
Parameter('Reactbn', 1)
Parameter('Reactbo', 1)
Parameter('Reactbp', 1)
Parameter('Reactbq', 1)

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
Rule('M1_binds_C1', M1(b=None, state='A') + C1(b=None, state='A') | M1(b=1, state='A') % C1(b=1, state='A'), Reactbc, Reactbg)
Rule('M2_binds_C1', M2(b=None, state='A') + C1(b=None, state='A') | M2(b=1, state='A') % C1(b=1, state='A'), Reactbd, Reactbh)
Rule('M8_binds_C1', M8(b=None, state='A') + C1(b=None, state='A') | M8(b=1, state='A') % C1(b=1, state='A'), Reactbe, Reactbi)
Rule('M9_binds_C1', M9(b=None, state='A') + C1(b=None, state='A') | M9(b=1, state='A') % C1(b=1, state='A'), Reactbf, Reactbj)

Rule('M1_deg_C1', M1(b=1) % C1(b=1, state='A') >> M1(b=None, state='A') + C1(b=None, state='D'), Reactaa)
Rule('M2_deg_C1', M2(b=1) % C1(b=1, state='A') >> M2(b=None, state='A') + C1(b=None, state='D'), Reactab)
Rule('M8_deg_C1', M8(b=1) % C1(b=1, state='A') >> M8(b=None, state='A') + C1(b=None, state='D'), Reactac)
Rule('M9_deg_C1', M9(b=1) % C1(b=1, state='A') >> M9(b=None, state='A') + C1(b=None, state='D'), Reactad)

Rule('M1_binds_T1', M1(b=None, state='A') + T1(b=None, state='A') | M1(b=1, state='A') % T1(b=1, state='A'), Reactae, Reactaf)
Rule('M2_binds_T1', M2(b=None, state='A') + T1(b=None, state='A') | M2(b=1, state='A') % T1(b=1, state='A'), Reactag, Reactah)
Rule('M8_binds_T1', M8(b=None, state='A') + T1(b=None, state='A') | M8(b=1, state='A') % T1(b=1, state='A'), Reactai, Reactaj)
Rule('M9_binds_T1', M9(b=None, state='A') + T1(b=None, state='A') | M9(b=1, state='A') % T1(b=1, state='A'), Reactak, Reactal)
Rule('M1_binds_T2', M1(b=None, state='A') + T2(b=None, state='A') | M1(b=1, state='A') % T2(b=1, state='A'), Reactam, Reactan)
Rule('M2_binds_T2', M2(b=None, state='A') + T2(b=None, state='A') | M2(b=1, state='A') % T2(b=1, state='A'), Reactao, Reactap)
Rule('M8_binds_T2', M8(b=None, state='A') + T2(b=None, state='A') | M8(b=1, state='A') % T2(b=1, state='A'), Reactaq, Reactar)
Rule('M9_binds_T2', M9(b=None, state='A') + T2(b=None, state='A') | M9(b=1, state='A') % T2(b=1, state='A'), Reactas, Reactat)
Rule('M1_binds_T4', M1(b=None, state='A') + T4(b=None, state='A') | M1(b=1, state='A') % T4(b=1, state='A'), Reactau, Reactav)
Rule('M2_binds_T4', M2(b=None, state='A') + T4(b=None, state='A') | M2(b=1, state='A') % T4(b=1, state='A'), Reactaw, Reactax)
Rule('M8_binds_T4', M8(b=None, state='A') + T4(b=None, state='A') | M8(b=1, state='A') % T4(b=1, state='A'), Reactay, Reactaz)
Rule('M9_binds_T4', M9(b=None, state='A') + T4(b=None, state='A') | M9(b=1, state='A') % T4(b=1, state='A'), Reactba, Reactbb)

Rule('M1_inhib', M1(b=None, state='A') >> M1(b=None, state='I') , Reactbk)
Rule('M2_inhib', M2(b=None, state='A') >> M2(b=None, state='I') , Reactbl)
Rule('M8_inhib', M8(b=None, state='A') >> M8(b=None, state='I') , Reactbm)
Rule('M9_inhib', M9(b=None, state='A') >> M9(b=None, state='I') , Reactbn)
Rule('T1_inhib', T1(b=None, state='A') >> T1(b=None, state='I') , Reactbo)
Rule('T2_inhib', T2(b=None, state='A') >> T2(b=None, state='I') , Reactbp)
Rule('T4_inhib', T4(b=None, state='A') >> T4(b=None, state='I') , Reactbq)

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
