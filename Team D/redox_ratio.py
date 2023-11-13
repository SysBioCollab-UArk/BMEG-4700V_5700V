from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()
print("ll")

#create monomers

#reaction 1 monomers
Monomer("NADH",['b'])
Monomer("Ox",['b'])

#reaction 2 monomers
Monomer("NAD",['b'])

#initial conditions
Parameter("NADH_init",100) #find initial condition values
Parameter("Ox_init",50)    #find initial condition values
Parameter("NAD_init", 60)  #find initial condition values

#initialize monomers
Initial(NADH(b=None),NADH_init)
Initial(Ox(b=None),Ox_init)
Initial(NAD(b=None),NAD_init)

#rates
Parameter("kp1",0.01) #find rate values
Parameter("km1",1)    #find rate values
Parameter("kp2",0.01) #find rate values
Parameter("km2", 1)   #find rate values

#rules
Rule("NADH_binds_Ox",NADH(b=None) + Ox(b=None) | NADH(b=1) % Ox(b=1), kp1, km1)
Rule("NADH_to_NAD", NADH(b=1) % Ox(b=1) | NAD(b=None) + Ox(b=None), kp2, km2)

#create observables to plot
Observable("NADH_free", NADH(b=None))
Observable("Ox_free", Ox(b=None))
Observable("NADH_Ox_bound", NADH(b=1) % Ox(b=1))
Observable("NAD_free", NAD(b=None))

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
