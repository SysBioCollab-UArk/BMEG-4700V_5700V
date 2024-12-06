from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers
Monomer('Complex_I', ['fmn'])
Monomer('FMN', ['cI', 'e_'], {'e_': ['_0', '_1', '_2']})
Monomer('NADH')
Monomer('NADplus')
Monomer('Hplus')
Monomer('e_', ['loc'], {'loc': ['mat', 'mem']})
Monomer('O2')
Monomer('O2_') # superoxide
Monomer('Q', ['e_'], {'e_': ['_0', '_1', '_2']})
Monomer('QH2')



# Initial conditions
Parameter('CI_0', 100)
Parameter('NADH_0', 1e4)
Parameter('Q_0', 100)
Parameter('O2_0', 100)

Initial(Complex_I(fmn=1) % FMN(cI=1, e_='_0'), CI_0)
Initial(NADH(), NADH_0)
Initial(Q(e_='_0'), Q_0)
Initial(O2(), O2_0)



# Rules
Parameter('k_NADH_ox', 1)
Parameter('k_superox_FMN', 1)
Parameter('k_FMN_Q_reduc', 1)
Parameter('k_Q_QH2', 1)

Rule('NADH_oxidation',
     NADH() + FMN(cI=ANY, e_='_0') >> NADplus() + Hplus() + FMN(cI=ANY, e_='_2'), k_NADH_ox)

Rule('superoxide_FMN_2_1',
     O2() + FMN(cI=ANY, e_='_2') >> O2_() + FMN(cI=ANY, e_='_1'), k_superox_FMN)

Rule('superoxide_FMN_1_0',
     O2() + FMN(cI=ANY, e_='_1') >> O2_() + FMN(cI=ANY, e_='_0'), k_superox_FMN)

Rule('FMN_Q_reduction_2_0',
     FMN(cI=ANY, e_='_2') + Q(e_='_0') >> FMN(cI=ANY, e_='_1') + Q(e_='_1'), k_FMN_Q_reduc)

Rule('FMN_Q_reduction_1_0',
     FMN(cI=ANY, e_='_1') + Q(e_='_0') >> FMN(cI=ANY, e_='_0') + Q(e_='_1'), k_FMN_Q_reduc)

Rule('FMN_Q_reduction_2_1',
     FMN(cI=ANY, e_='_2') + Q(e_='_1') >> FMN(cI=ANY, e_='_1') + Q(e_='_2'), k_FMN_Q_reduc)

Rule('FMN_Q_reduction_1_1',
     FMN(cI=ANY, e_='_1') + Q(e_='_1') >> FMN(cI=ANY, e_='_0') + Q(e_='_2'), k_FMN_Q_reduc)

Rule('Q_reduction_QH2',
     Q(e_='_2') + Hplus() + Hplus() >> QH2(), k_Q_QH2)




# Observables
Observable('FMN_e0', FMN(e_='_0'))
Observable('FMN_e1', FMN(e_='_1'))
Observable('FMN_e2', FMN(e_='_2'))
Observable('Q_e0', Q(e_='_0'))
Observable('Q_e1', Q(e_='_1'))
Observable('Q_e2', Q(e_='_2'))
Observable('QH2_tot', QH2())
Observable('O2_tot', O2())
Observable('SuperOx', O2_())

# Simulation commands + plotting
tspan= np.linspace(0, 0.1, 101)
sim=ScipyOdeSimulator(model, tspan, verbose=True)
output=sim.run()

for obs in [o for o in model.observables if o.name not in ['O2_tot', 'SuperOx']]:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc=(0.5,0.55), ncol=2)
plt.tight_layout()

plt.figure()
for obs in [o for o in model.observables if o.name in ['O2_tot', 'SuperOx']]:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc='best')
plt.tight_layout()

plt.show()
