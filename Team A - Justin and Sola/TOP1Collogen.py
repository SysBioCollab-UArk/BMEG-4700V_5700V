from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt



Model()



Monomer('M1', ['b'])
Monomer('M2', ['b'])
Monomer('M8', ['b'])
Monomer('M9', ['b'])
Monomer('T1', ['b'])
Monomer('T2', ['b'])
Monomer('T4', ['b'])

Monomer('C1', ['b', 'state'], {'state': ['A', 'D']})

Initial(M1(b=None), Parameter('M1', 750))
Initial(M2(b=None), Parameter('M2', 750))
Initial(M8(b=None), Parameter('M8', 750))
Initial(M9(b=None), Parameter('M9', 750))
Initial(T1(b=None), Parameter('T1', 750))
Initial(T2(b=None), Parameter('T2', 750))
Initial(T4(b=None), Parameter('T4', 750))

Initial(C1(b=None, state='A'), Parameter('C1', 1500))



# collegen degridation
Parameter('Reactaa', 500)
Parameter('Reactab', 500)
Parameter('Reactac', 500)
Parameter('Reactad', 500)

# Forward parameters
Parameter('Reactae', 500)
Parameter('Reactag', 500)
Parameter('Reactai', 500)
Parameter('Reactak', 500)
Parameter('Reactam', 500)
Parameter('Reactao', 500)
Parameter('Reactaq', 500)
Parameter('Reactas', 500)
Parameter('Reactau', 500)
Parameter('Reactaw', 500)
Parameter('Reactay', 500)
Parameter('Reactba', 500)
# Reverse parameters
Parameter('Reactaf', 500)
Parameter('Reactah', 500)
Parameter('Reactaj', 500)
Parameter('Reactal', 500)
Parameter('Reactan', 500)
Parameter('Reactap', 500)
Parameter('Reactar', 500)
Parameter('Reactat', 500)
Parameter('Reactav', 500)
Parameter('Reactax', 500)
Parameter('Reactaz', 500)
Parameter('Reactbb', 500)










Rule('M1_deg_C1', M1(b=None) + C1(b=None, state='A') >> M1(b=None) + C1(b=None, state='D'), Reactaa)
Rule('M2_deg_C1', M2(b=None) + C1(b=None, state='A') >> M2(b=None) + C1(b=None, state='D'), Reactab)
Rule('M8_deg_C1', M8(b=None) + C1(b=None, state='A') >> M8(b=None) + C1(b=None, state='D'), Reactac)
Rule('M9_deg_C1', M9(b=None) + C1(b=None, state='A') >> M9(b=None) + C1(b=None, state='D'), Reactad)

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


Observable('unbound_M1', M1(b=None))
Observable('unbound_T1', T1(b=None))

Observable('bound_M1', M1(b=1))
Observable('bound_T1', T1(b=1))

Observable('decC1', Lacto(b=None, state='D'))
Observable('actC1', Lacto(b=None, state='A'))



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














