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
Initial(M1(b=None), Parameter('M1_init', 2.35))
Initial(M2(b=None), Parameter('M2_init', 363.12))
Initial(M8(b=None), Parameter('M8_init', 2.08))
Initial(M9(b=None), Parameter('M9_init', 85.57))
Initial(T1(b=None), Parameter('T1_init', 75.16))
Initial(T2(b=None), Parameter('T2_init', 81.74))
Initial(T4(b=None), Parameter('T4_init', 0.97))

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
# Reverse parameters
#Parameter('Reactaf', .0001)
#Parameter('Reactah', .0001)
#Parameter('Reactaj', .0001)
#Parameter('Reactal', .0001)
#Parameter('Reactan', .0001)
#Parameter('Reactap', .0001)
#Parameter('Reactar', .0001)
#Parameter('Reactat', .0001)
#Parameter('Reactav', .0001)
#Parameter('Reactax', .0001)
#Parameter('Reactaz', .0001)
#Parameter('Reactbb', .0001)

# Rules
Rule('M1_deg_C1', M1(b=None) + C1(b=None, state='A') >> M1(b=None) + C1(b=None, state='D'), Reactaa)
Rule('M2_deg_C1', M2(b=None) + C1(b=None, state='A') >> M2(b=None) + C1(b=None, state='D'), Reactab)
Rule('M8_deg_C1', M8(b=None) + C1(b=None, state='A') >> M8(b=None) + C1(b=None, state='D'), Reactac)
Rule('M9_deg_C1', M9(b=None) + C1(b=None, state='A') >> M9(b=None) + C1(b=None, state='D'), Reactad)

Rule('M1_binds_T1', M1(b=None) + T1(b=None) >> M1(b=1) % T1(b=1), Reactae)
Rule('M2_binds_T1', M2(b=None) + T1(b=None) >> M2(b=1) % T1(b=1), Reactag)
Rule('M8_binds_T1', M8(b=None) + T1(b=None) >> M8(b=1) % T1(b=1), Reactai)
Rule('M9_binds_T1', M9(b=None) + T1(b=None) >> M9(b=1) % T1(b=1), Reactak)
Rule('M1_binds_T2', M1(b=None) + T2(b=None) >> M1(b=1) % T2(b=1), Reactam)
Rule('M2_binds_T2', M2(b=None) + T2(b=None) >> M2(b=1) % T2(b=1), Reactao)
Rule('M8_binds_T2', M8(b=None) + T2(b=None) >> M8(b=1) % T2(b=1), Reactaq)
Rule('M9_binds_T2', M9(b=None) + T2(b=None) >> M9(b=1) % T2(b=1), Reactas)
Rule('M1_binds_T4', M1(b=None) + T4(b=None) >> M1(b=1) % T4(b=1), Reactau)
Rule('M2_binds_T4', M2(b=None) + T4(b=None) >> M2(b=1) % T4(b=1), Reactaw)
Rule('M8_binds_T4', M8(b=None) + T4(b=None) >> M8(b=1) % T4(b=1), Reactay)
Rule('M9_binds_T4', M9(b=None) + T4(b=None) >> M9(b=1) % T4(b=1), Reactba)

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



#Observable('bound_M1', M1(b=ANY))
#Observable('bound_M2', M2(b=ANY))
#Observable('bound_M8', M8(b=ANY))
#Observable('bound_M9', M9(b=ANY))
#
#
#Observable('bound_T1', T1(b=ANY))
#Observable('bound_T2', T2(b=ANY))
#Observable('bound_T4', T4(b=ANY))



# Simulation commands
tspan = np.linspace(0, 3, 501)
sim = ScipyOdeSimulator(model, tspan, verbose=True)

print(len(model.species))
for i, sp in enumerate (model.species):
    print(i,sp)
print()
for i, ode in enumerate(model.odes):
    print(model.species[i],ode)




#quit()


output = sim.run()

for obs in model.observables:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc=0)

plt.show()

