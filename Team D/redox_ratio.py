from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()
print("ll")

#create monomers
Monomer("NADH",['b'])
Monomer("Ox",['b'])

#initial conditions
Parameter("NADH_init",100)
Parameter("Ox_init",50)
Initial(NADH(b=None),NADH_init)
Initial(Ox(b=None),Ox_init)

#rates
Parameter("kp1",1)
Parameter("km1",10)

#rules
Rule("NADH_binds_Ox",NADH(b=None) + Ox(b=None) | NADH(b=1) % Ox(b=1), kp1, km1)

#create observables to plot
Observable("NADH_free", NADH(b=None))
Observable("Ox_free", Ox(b=None))
Observable("NADH_Ox_bound", NADH(b=1) % Ox(b=1))

print(model)
print(model.monomers)
print(model.observables)
#quit()

## simulation commands ##
tspan = np.linspace(0,1,101)
sim = ScipyOdeSimulator(model, tspan, verbose = True)

output = sim.run()


## plotting
for obs in model.observables:
    plt.plot(tspan,output.observables[obs.name],lw = 2, label = obs.name)
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend(loc = 0)
plt.show()
