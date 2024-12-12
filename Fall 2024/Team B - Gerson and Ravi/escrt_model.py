from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers
Monomer('cargo', ['escrt0', 'state'], {'state': ['x', 'ub']})
Monomer('ESCRT_0', ['cargo', 'escrt1'])
Monomer('ESCRT_I', ['escrt0', 'escrt2'])
Monomer('ESCRT_II', ['escrt1', 'escrt2', 'escrt3'])
Monomer('ESCRT_III', ['escrt2', 'escrt2', 'vesicle'])
Monomer('vesicle', ['escrt3', 'cargo', 'cargo'])
Monomer('VPS4_VTA1', ['escrt3', 'escrt3'])

# Initial conditions
Parameter('cargo_init', 1e5)
Parameter('ESCRT_0_init', 100)
Parameter('ESCRT_I_init', 100)
Parameter('ESCRT_II_init', 100)
Parameter('ESCRT_III_init', 100)
Parameter('VPS4_VTA1_init', 100)

Initial(cargo(escrt0=None, state='x'), cargo_init)
Initial(ESCRT_0(cargo=None, escrt1=None), ESCRT_0_init)
Initial(ESCRT_I(escrt0=None, escrt2=None), ESCRT_I_init)
Initial(ESCRT_II(escrt1=None, escrt2=None, escrt3=None), ESCRT_II_init)
Initial(ESCRT_III(escrt2=MultiState(None, None), vesicle=None), ESCRT_III_init)
Initial(VPS4_VTA1(escrt3=MultiState(None, None)), VPS4_VTA1_init)

# Rules
Parameter('k0', 1)
Rule('cargo_ubiquitinated',
     cargo(escrt0=None, state='x') >> cargo(escrt0=None, state='ub'), k0)

Parameter('kf1', 1)
Parameter('kr1', 1)
Rule('ESCRT0_binds_ubiq_cargo',
     ESCRT_0(cargo=None, escrt1=None) + cargo(escrt0=None, state='ub') |
     ESCRT_0(cargo=1, escrt1=None) % cargo(escrt0=1, state='ub'), kf1, kr1)

Parameter('kf2', 1)
Parameter('kr2', 1)
Rule('ESCRT0_recruits_ESCRT1',
     ESCRT_0(cargo=1, escrt1=None) % cargo(escrt0=1, state='ub') + ESCRT_I(escrt0=None, escrt2=None) |
     ESCRT_0(cargo=1, escrt1=2) % cargo(escrt0=1, state='ub') % ESCRT_I(escrt0=2, escrt2=None), kf2, kr2)

Parameter('kf3', 1)
Parameter('kr3', 1)
Rule('ESCRT1_recruits_ESCRT2',
     ESCRT_0(cargo=1, escrt1=2) % cargo(escrt0=1, state='ub') % ESCRT_I(escrt0=2, escrt2=None) +
     ESCRT_II(escrt1=None, escrt2=None, escrt3=None) |
     ESCRT_0(cargo=1, escrt1=2) % cargo(escrt0=1, state='ub') % ESCRT_I(escrt0=2, escrt2=3) %
     ESCRT_II(escrt1=3, escrt2=None, escrt3=None), kf3, kr3)

Parameter('kf4', 1)
Parameter('kr4', 1)
Rule('Membrane_invagination',
     ESCRT_II(escrt1=ANY, escrt2=None, escrt3=None) + ESCRT_II(escrt1=ANY, escrt2=None, escrt3=None) |
     ESCRT_II(escrt1=ANY, escrt2=1, escrt3=None) % ESCRT_II(escrt1=ANY, escrt2=1, escrt3=None), kf4, kr4)

Parameter('k5', 0.1)
Rule('Vesicle_formation',
     ESCRT_II(escrt1=ANY, escrt2=1, escrt3=None) % ESCRT_II(escrt1=ANY, escrt2=1, escrt3=None) +
     ESCRT_III(escrt2=MultiState(None, None), vesicle=None) >>
     ESCRT_II(escrt1=ANY, escrt2=1, escrt3=2) % ESCRT_II(escrt1=ANY, escrt2=1, escrt3=3) %
     ESCRT_III(escrt2=MultiState(2, 3), vesicle=4) % vesicle(escrt3=4, cargo=MultiState(None, None)), k5)

Parameter('k6', 1)
Rule('Cargo_deubiquitination',
     ESCRT_III(escrt2=MultiState(ANY, ANY), vesicle=1) % vesicle(escrt3=1, cargo=None) %
     cargo(escrt0=2, state='ub') % ESCRT_0(cargo=2, escrt1=ANY) >>
     ESCRT_III(escrt2=MultiState(ANY, ANY), vesicle=1) % vesicle(escrt3=1, cargo=2) %
     cargo(escrt0=2, state='x') % ESCRT_0(cargo=None, escrt1=ANY), k6)

Parameter('k7', 1)
Rule('ESCRT_0_release',
     ESCRT_III(escrt2=MultiState(ANY, ANY)) % ESCRT_I(escrt0=1, escrt2=ANY) % ESCRT_0(cargo=None, escrt1=1) >>
     ESCRT_III(escrt2=MultiState(ANY, ANY)) % ESCRT_I(escrt0=None, escrt2=ANY) + ESCRT_0(cargo=None, escrt1=None),
     k7)

Parameter('k8', 100)
Rule('VPS4_VTA1_binds_ESCRT3',
     VPS4_VTA1(escrt3=MultiState(None, None)) +
     ESCRT_III(escrt2=MultiState(2, 3), vesicle=ANY) %
     ESCRT_I(escrt0=None, escrt2=4) % ESCRT_II(escrt1=4, escrt2=1, escrt3=2) %
     ESCRT_II(escrt1=5, escrt2=1, escrt3=3) % ESCRT_I(escrt0=None, escrt2=5)
     >>
     VPS4_VTA1(escrt3=MultiState(1, 2)) %
     ESCRT_III(escrt2=MultiState(1, 2), vesicle=ANY) +
     ESCRT_I(escrt0=None, escrt2=None) + ESCRT_II(escrt1=None, escrt2=None, escrt3=None) +
     ESCRT_II(escrt1=None, escrt2=None, escrt3=None) + ESCRT_I(escrt0=None, escrt2=None),
     k8)

Parameter('k9', 1)
Rule('Vesicle_release',
     VPS4_VTA1(escrt3=MultiState(1, 2)) % ESCRT_III(escrt2=MultiState(1, 2), vesicle=3) %
     vesicle(escrt3=3) >>
     VPS4_VTA1(escrt3=MultiState(None, None)) + vesicle(escrt3=None), k9)

# Observables
# Observable('free_ESCRT0', ESCRT_0(cargo=None, escrt1=None))
Observable('cargo_ESCRT0', cargo(escrt0=1) % ESCRT_0(cargo=1, escrt1=None))
Observable('cargo_ESCRT0_ESCRT1', cargo(escrt0=1) % ESCRT_0(cargo=1, escrt1=ANY))
Observable('cargo_ESCRT0_ESCRT1_ESCRT2',
           cargo(escrt0=1) % ESCRT_0(cargo=1, escrt1=2) % ESCRT_I(escrt0=2, escrt2=ANY))
Observable('initial_bud', ESCRT_II(escrt2=ANY) % cargo(state='ub'))
Observable('ESCRT3_assembly_cargo_ub',
           ESCRT_III(escrt2=MultiState(ANY, ANY)) % cargo(escrt0=ANY, state='ub'))
Observable('ESCRT3_assembly_cargo_dub',
           ESCRT_III(escrt2=MultiState(ANY, ANY)) % cargo(escrt0=ANY, state='x'))
Observable('neck_constriction', VPS4_VTA1(escrt3=MultiState(ANY, ANY)))
Observable('total_vesicles', vesicle())
Observable('free_vesicles', vesicle(escrt3=None))

# simulation commands
tspan = np.linspace(0, 10, 1001)
sim = ScipyOdeSimulator(model, tspan, verbose=True)

# for i, sp in enumerate(model.species):
#     print(i, sp)

output = sim.run()

# Plots of species at different steps along the ESCRT pathway
plt.figure(figsize=(6.4*1.5, 4.8*1.2))
for obs in model.observables:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()

# Plots of ubiquitinated and deubiquitinated cargo before and after vesicle formation and number of free vesicles
plt.figure()
plt.plot(tspan, output.observable(cargo(escrt0=1, state='ub') % ESCRT_0(cargo=1)), lw=2,
         label='cargo_ub_bound_ESCRT')
plt.plot(tspan, output.observable(cargo(escrt0=1, state='x') % vesicle(cargo=1)), lw=2,
         label='cargo_x_bound_vesicle')
plt.plot(tspan, output.observables['free_vesicles'], lw=2, label='free_vesicles')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc='best')
plt.tight_layout()

# Plot of ESCRT III concentration (goes to zero over time)
plt.figure()
plt.plot(tspan, output.observable(ESCRT_III()), lw=2, label='ESCRT_III')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc='best')
plt.tight_layout()

# Plots of total, deubiquitinated, and ubiquitinated cargo
plt.figure()
plt.plot(tspan, output.observable(cargo()), lw=2, label='cargo_total')
plt.plot(tspan, output.observable(cargo(state='x')), lw=2, label='cargo_dUb')
plt.plot(tspan, output.observable(cargo(state='ub')), lw=2, label='cargo_ub')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc='best')
plt.tight_layout()

plt.show()
