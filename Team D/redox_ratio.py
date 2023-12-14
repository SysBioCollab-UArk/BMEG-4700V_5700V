from pysb import Model, Monomer, Parameter, Initial, Rule, Observable, Expression
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

# Create a new model instance
model = Model("RedoxModel")

# Define Monomers
Monomer("NADH", ['state', 'b'], {'state': ['free', 'bound']})  # Added 'b' site
Monomer("Ox", ['b'])  # Oxidase
Monomer("NAD_plus", ['b'])
Monomer("Inhibitor", ['type'], {'type': ['oxamate', 'rot', 'oligo', 'fccp']})

# Define Parameters
Parameter("NADH_free_init", 100)
Parameter("NADH_bound_init", 50)
Parameter("Ox_init", 50)
Parameter("NAD_plus_init", 0)
Parameter("Inhibitor_init", 0)  # Initial amount of inhibitor
Parameter("alpha", 0.01)  # Effective rate of NADH oxidation
Parameter("beta_eq", 1)   # Effective rate of NAD+ reduction
Parameter("inhibition_factor", 0.5)  # Effect of inhibitors on the reaction rates

# Define Expressions for rates
Expression("alpha_inhibited", alpha * inhibition_factor)
Expression("beta_eq_inhibited", beta_eq * inhibition_factor)

# Define Initials
Initial(NADH(state='free', b=None), NADH_free_init)
Initial(NADH(state='bound', b=None), NADH_bound_init)
Initial(Ox(b=None), Ox_init)
Initial(NAD_plus(b=None), NAD_plus_init)
Initial(Inhibitor(type='oxamate'), Inhibitor_init)
Initial(Inhibitor(type='rot'), Inhibitor_init)
Initial(Inhibitor(type='oligo'), Inhibitor_init)
Initial(Inhibitor(type='fccp'), Inhibitor_init)

# Define Reaction Rules
Rule('NADH_binding',
     NADH(state='free', b=None) + Ox(b=None) >>
     NADH(state='bound', b=1) % Ox(b=1),
     Parameter('k_bind', 1e-3))

Rule('NADH_unbinding',
     NADH(state='bound', b=1) % Ox(b=1) >>
     NADH(state='free', b=None) + Ox(b=None),
     Parameter('k_unbind', 1e-3))

for inhibitor_type in ['oxamate', 'rot', 'oligo', 'fccp']:
    Rule(f'Inhibit_NADH_free_{inhibitor_type}',
         NADH(state='free') + Inhibitor(type=inhibitor_type) >>
         NADH(state='free') % Inhibitor(type=inhibitor_type),
         Parameter(f'k_inhibit_free_{inhibitor_type}', 1e-4))

    Rule(f'Inhibit_NADH_bound_{inhibitor_type}',
         NADH(state='bound') + Inhibitor(type=inhibitor_type) >>
         NADH(state='bound') % Inhibitor(type=inhibitor_type),
         Parameter(f'k_inhibit_bound_{inhibitor_type}', 1e-4))

# Define Observables
Observable("NADH_free", NADH(state='free'))
Observable("NADH_bound", NADH(state='bound'))

# Simulation setup
tspan = np.linspace(0, 100, 101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
output = sim.run()

# Calculate NADH ratio
nadh_ratio = output.observables['NADH_free'] / output.observables['NADH_bound']

# Plotting results
plt.figure(figsize=(10, 6))
plt.plot(tspan, output.observables['NADH_free'], label='NADH_free')
plt.plot(tspan, output.observables['NADH_bound'], label='NADH_bound')
plt.plot(tspan, nadh_ratio, label='NADH_ratio', linestyle='--')
plt.xlabel("Time")
plt.ylabel("Concentration or Ratio")
plt.legend(loc='best')
plt.show()
