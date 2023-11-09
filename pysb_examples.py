from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

Monomer('P')
Observable('P_tot', P())
Parameter('k1', 1)
Parameter('k2', 50)
Expression('rate_P_exp', k1/(k2+P_tot))
Rule('Protein_expression', None >> P(), rate_P_exp)

Monomer('E')
Parameter('E_0', 100)
Initial(E(), E_0)
Monomer('F')
Observable('F_tot', F())
Parameter('k3', 1)
Parameter('k4', 100)
Expression('rate_fragment', k3/(k4+F_tot))
Rule('Protein_fragmentation', E() + P() >> E() + F() + F(), rate_fragment)

Parameter('k_F_deg', 0.1)
Rule('Fragment_degradation', F() >> None, k_F_deg)

tspan = np.linspace(0, 100, 1001)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
output = sim.run()

for obs in model.observables:
    plt.plot(tspan, output.observables[obs.name], label=obs.name)
plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc=0)
plt.tight_layout()

plt.show()
