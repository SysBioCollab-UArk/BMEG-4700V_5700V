from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()
#pysb automatically creates variables behind the scenes when you call these functions
# rules = interaction, parameters - rates of rxn, initial amts of species, expressions - rxns that occur not with constant rates but rates that vary depending on concentration of molecules, more complex kinetics. Compartments - cells e.g. extra vs intracellular, nuclear


# Monomers
Monomer('A',['b'])
Monomer('B',['a','state'],{'state':['u','p']})

#print(A)
#print(type(A))

# Initial conditions
Parameter('A_0',100)
Parameter('B_0',80)
Initial(A(b=None), A_0)
Initial(B(a=None,state='u'), B_0)

# Rules
Parameter('kf_AB',1)
Parameter('kr_AB',10)
Rule('A_B_binds', A(b=None) + B(a=None,state='u') | A(b=1) % B(a=1,state='u'), kf_AB, kr_AB)
# Reaction center - thing that changes in the rule
# Reaction context - doesnt change, but required condition for rule to occur
# "|" means reversible
# '%' means that they are in the same complex with each other, not that they are bound to each other
# 'ring closure' vs free floating model. '+' means they are not in the same complex
# Mechanistic model vs statistical (ML) model

Parameter('k_B_phos',0.1)
Rule('AB_phosphorylate', A(b=1) % B(a=1,state='u') >> A(b=1) % B(a=1,state='p'), k_B_phos)

Parameter('k_ABp_unbind',1)
Rule('A_Bp_unbinds', A(b=1) % B(a=1,state='p') >> A(b=None) + B(a=None,state='p'), k_ABp_unbind)

Parameter('k_Bp_unphos',10)
Rule('Bp_unphos', B(a=None, state='p') >> B(a=None, state='u'), k_Bp_unphos)

# Observables
# patterns; match species
Observable('B_u', B(state='u'))
Observable('B_p', B(state='p'))
Observable('Free_A', A(b=None))
Observable('Free_B', B(a=None))
Observable('AB_complex', A(b=1) % B(a=1))
Observable('ABp_complex', A(b=1) % B(a=1, state='p'))


# Simulation commands + plotting
tspan = np.linspace(0,0.1,101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
result = sim.run()

for obs in model.observables:
    plt.plot(tspan, result.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc=0, ncol=2)
plt.tight_layout()

plt.show()

