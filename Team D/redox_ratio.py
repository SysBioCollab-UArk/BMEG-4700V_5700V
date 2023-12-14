from pysb import Model, Monomer, Parameter, Initial, Rule, Observable
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

model = Model("RedoxModel")

Monomer('NADH', ['state', 'b'], {'state': ['free', 'bound', 'inhibited']})
Monomer('Ox', ['b'])
Monomer('NAD_plus', ['b'])
Monomer('Inhibitor', ['type', 'b'], {'type': ['oxamate', 'rot', 'oligo', 'fccp']})

Parameter("NADH_free_init", 100)
Parameter("NADH_bound_init", 50)
Parameter("Ox_init", 50)
Parameter("NAD_plus_init", 0)
Parameter("Inhibitor_init", 0)
Parameter("k_inhibit_free_oxamate", 1e-4)
Parameter("k_inhibit_bound_oxamate", 1e-4)
Parameter("k_reverse_inhibit_free_oxamate", 1e-4)
Parameter("k_inhibit_free_rot", 1e-4)
Parameter("k_inhibit_bound_rot", 1e-4)
Parameter("k_reverse_inhibit_free_rot", 1e-4)
Parameter("k_inhibit_free_oligo", 1e-4)
Parameter("k_inhibit_bound_oligo", 1e-4)
Parameter("k_reverse_inhibit_free_oligo", 1e-4)
Parameter("k_inhibit_free_fccp", 1e-4)
Parameter("k_inhibit_bound_fccp", 1e-4)
Parameter("k_reverse_inhibit_free_fccp", 1e-4)
Parameter("k_bind", 1e-3)
Parameter("k_unbind", 1e-3)
Parameter("alpha", 0.01)
Parameter("beta_eq", 1)
Parameter("inhibition_factor", 0.5)

Initial(NADH(state='free', b=None), NADH_free_init)
Initial(NADH(state='bound', b=None), NADH_bound_init)
Initial(Ox(b=None), Ox_init)
Initial(NAD_plus(b=None), NAD_plus_init)
Initial(Inhibitor(type='oxamate', b=None), Inhibitor_init)
Initial(Inhibitor(type='rot', b=None), Inhibitor_init)
Initial(Inhibitor(type='oligo', b=None), Inhibitor_init)
Initial(Inhibitor(type='fccp', b=None), Inhibitor_init)

# Binding/unbinding rules
Rule('NADH_binding',
     NADH(state='free', b=None) + Ox(b=None) >>
     NADH(state='bound', b=1) % Ox(b=1),
     k_bind)

Rule('NADH_unbinding',
     NADH(state='bound', b=1) % Ox(b=1) >>
     NADH(state='free', b=None) + Ox(b=None),
     k_unbind)

# oxamate inhibition

Rule('Inhibit_NADH_free_oxamate',
     NADH(state='free', b=None) + Inhibitor(type='oxamate', b=None) >>
     NADH(state='inhibited', b=1) % Inhibitor(type='oxamate', b=1),
     k_inhibit_free_oxamate)

Rule('Inhibit_NADH_bound_oxamate',
     NADH(state='bound', b=1) % Ox(b=1) + Inhibitor(type='oxamate', b=None) >>
     NADH(state='inhibited', b=2) % Inhibitor(type='oxamate', b=2),
     k_inhibit_bound_oxamate)

Rule('Reverse_Inhibit_NADH_free_oxamate',
     NADH(state='inhibited', b=1) % Inhibitor(type='oxamate', b=1) >>
     NADH(state='free', b=None) + Inhibitor(type='oxamate', b=None),
     k_reverse_inhibit_free_oxamate)

# Rot inhibition
Rule('Inhibit_NADH_free_rot',
     NADH(state='free', b=None) + Inhibitor(type='rot', b=None) >>
     NADH(state='inhibited', b=1) % Inhibitor(type='rot', b=1),
     k_inhibit_free_rot)

Rule('Inhibit_NADH_bound_rot',
     NADH(state='bound', b=1) % Ox(b=1) + Inhibitor(type='rot', b=None) >>
     NADH(state='inhibited', b=2) % Inhibitor(type='rot', b=2),
     k_inhibit_bound_rot)

Rule('Reverse_Inhibit_NADH_free_rot',
     NADH(state='inhibited', b=1) % Inhibitor(type='rot', b=1) >>
     NADH(state='free', b=None) + Inhibitor(type='rot', b=None),
     k_reverse_inhibit_free_rot)

# Oligomycin inhibition
Rule('Inhibit_NADH_free_oligo',
     NADH(state='free', b=None) + Inhibitor(type='oligo', b=None) >>
     NADH(state='inhibited', b=1) % Inhibitor(type='oligo', b=1),
     k_inhibit_free_oligo)

Rule('Inhibit_NADH_bound_oligo',
     NADH(state='bound', b=1) % Ox(b=1) + Inhibitor(type='oligo', b=None) >>
     NADH(state='inhibited', b=2) % Inhibitor(type='oligo', b=2),
     k_inhibit_bound_oligo)

Rule('Reverse_Inhibit_NADH_free_oligo',
     NADH(state='inhibited', b=1) % Inhibitor(type='oligo', b=1) >>
     NADH(state='free', b=None) + Inhibitor(type='oligo', b=None),
     k_reverse_inhibit_free_oligo)

# FCCP inhibition
Rule('Inhibit_NADH_free_fccp',
     NADH(state='free', b=None) + Inhibitor(type='fccp', b=None) >>
     NADH(state='inhibited', b=1) % Inhibitor(type='fccp', b=1),
     k_inhibit_free_fccp)

Rule('Inhibit_NADH_bound_fccp',
     NADH(state='bound', b=1) % Ox(b=1) + Inhibitor(type='fccp', b=None) >>
     NADH(state='inhibited', b=2) % Inhibitor(type='fccp', b=2),
     k_inhibit_bound_fccp)

Rule('Reverse_Inhibit_NADH_free_fccp',
     NADH(state='inhibited', b=1) % Inhibitor(type='fccp', b=1) >>
     NADH(state='free', b=None) + Inhibitor(type='fccp', b=None),
     k_reverse_inhibit_free_fccp)

Observable("NADH_free", NADH(state='free'))
Observable("NADH_bound", NADH(state='bound'))

tspan = np.linspace(0, 100, 101)
sim = ScipyOdeSimulator(model, tspan)
results = sim.run()

# Calculate NADH ratio
nadh_ratio = results.observables['NADH_free'] / results.observables['NADH_bound']

# Plotting results
plt.figure(figsize=(10, 6))
plt.plot(tspan, results.observables['NADH_free'], label='NADH_free')
plt.plot(tspan, results.observables['NADH_bound'], label='NADH_bound')
plt.plot(tspan, nadh_ratio, label='NADH_ratio', linestyle='--')
plt.xlabel("Time")
plt.ylabel("Concentration or Ratio")
plt.legend(loc='best')
plt.show()
