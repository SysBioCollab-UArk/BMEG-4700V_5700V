from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()
# Monomers

Monomer('HIF1', ['p1', 'p2', 'a1'],
        {'p1': ['u', 'oh', 'ub'], 'p2': ['u', 'oh', 'ub'], 'a1': ['u', 'oh']})
Monomer('HIF2', ['p1','p2', 'a1'],
        {'p1': ['u', 'oh', 'ub'], 'p2': ['u', 'oh', 'ub'], 'a1': ['u', 'oh']})
Monomer('PHD1', ['hif_p'])
Monomer('PHD2', ['hif_p'])
Monomer('PHD3', ['hif_p'])
Monomer('FIH', ['hif_a1'])
Monomer('VHL', ['hif_p'])
Monomer('p300', ['hif_a1'])
Monomer('proteosome', ['hif_p'])

# Parameters
Parameter('kf_PHD3_binds_HIF1_p2')
Parameter('kr_PHD3_binds_HIF1_p2')
Parameter('k_PHD3_HIF1_hydroxy_p2')

Parameter('kf_PHD3_binds_HIF1_p1')
Parameter('kr_PHD3_binds_HIF1_p1')
Parameter('k_PHD3_HIF1_hydroxy_p1')

Parameter('kf_PHD2_binds_HIF1_p2')
Parameter('kr_PHD2_binds_HIF1_p2')
Parameter('k_PHD2_HIF1_hydroxy_p2')

Parameter('kf_PHD2_binds_HIF1_p1')
Parameter('kr_PHD2_binds_HIF1_p1')
Parameter('k_PHD2_HIF1_hydroxy_p1')

Parameter('kf_VHL_binds_HIF_1_p2')
Parameter('kr_VHL_binds_HIF_1_p2')
Parameter('k_VHL_ubiq_HIF_1_p2')

Parameter('kf_VHL_binds_HIF_1_p1')
Parameter('kr_VHL_binds_HIF_1_p1')
Parameter('k_VHL_ubiq_HIF_1_p1')

Parameter('kf_HIF_1_p1_proteo_binding')
Parameter('kr_HIF_1_p1_proteo_binding')
Parameter('k_HIF_1_p1_degrades')

Parameter('kf_HIF_1_p2_proteo_binding')
Parameter('kr_HIF_1_p2_proteo_binding')
Parameter('k_HIF_1_p2_degrades')

Parameter('kf_HIF1_binds_FIH')
Parameter('kr_HIF1_binds_FIH')
Parameter('k_HIF_1_hydroxylated_at_a1')

Parameter('kf_HIF1_binds_p300')
Parameter('kr_HIF1_binds_p300')

# Rules

#PHD3 hydroxylates HIF-1a at Proline 564
Rule('PHD3_binds_HIF1_p2', HIF1(p2='u') + PHD3(hif_p=None) | HIF1(p2=('u', 1)) % PHD3(hif_p=1),
     kf_PHD3_binds_HIF1_p2, kr_PHD3_binds_HIF1_p2 )
Rule('PHD3_hydroxy_HIF1_p2', HIF1(p2=('u', 1)) % PHD3(hif_p=1) >> HIF1(p2='oh') + PHD3(hif_p=None),
     k_PHD3_HIF1_hydroxy_p2)

#PHD3 hydroxylated HIF-1a at Proline 402
Rule('PHD3_binds_HIF1_p1', HIF1(p1='u') + PHD3(hif_p=None) | HIF1(p1=('u', 1)) % PHD3(hif_p=1),
     kf_PHD3_binds_HIF1_p1, kr_PHD3_binds_HIF1_p1)
Rule('PHD3_hydroxy_HIF1_p1', HIF1(p1=('u', 1)) % PHD3(hif_p=1) >> HIF1(p1='oh') + PHD3(hif_p=None),
     k_PHD3_HIF1_hydroxy_p1)

#PHD2 hydroxylates HIF-1a at Proline 564
Rule('PHD2_binds_HIF1_p2', HIF1(p2='u') + PHD2(hif_p=None) | HIF1(p2=('u', 1)) % PHD2(hif_p=1),
     kf_PHD2_binds_HIF1_p2, kr_PHD2_binds_HIF1_p2 )
Rule('PHD2_hydroxy_HIF1_p2', HIF1(p2=('u', 1)) % PHD2(hif_p=1) >> HIF1(p2='oh') + PHD2(hif_p=None),
     k_PHD2_HIF1_hydroxy_p2)

#PHD2 hydroxylates HIF-1a at Proline 402
Rule('PHD2_binds_HIF1_p1', HIF1(p1='u') + PHD2(hif_p=None) | HIF1(p1=('u', 1)) % PHD2(hif_p=1),
     kf_PHD2_binds_HIF1_p1, kr_PHD2_binds_HIF1_p1)
Rule('PHD2_hydroxy_HIF1_p1', HIF1(p1=('u', 1)) % PHD2(hif_p=1) >> HIF1(p1='oh') + PHD2(hif_p=None),
     k_PHD2_HIF1_hydroxy_p1)

#VHL binds to HIF-1 with 1-OH on p2 and ubiquinates
Rule('VHL_binds_HIF_1_p2', HIF1(p2='oh') + VHL(hif_p=None) | HIF1(p2=('oh', 1)) % VHL(hif_p=1),
     kf_VHL_binds_HIF_1_p2, kr_VHL_binds_HIF_1_p2)
Rule('VHL_ubiqu_HIF_1_p2',HIF1(p2=('oh', 1)) % VHL(hif_p=1) >> HIF1(p2 = 'ub') + VHL(hif_p=None),
     k_VHL_ubiq_HIF_1_p2)

#VHL binds to HIF-1 with 1-OH on p1 and ubiquinates
Rule('VHL_binds_HIF_1_p1', HIF1(p1='oh') + VHL(hif_p=None) | HIF1(p1=('oh', 1)) % VHL(hif_p=1),
     kf_VHL_binds_HIF_1_p1, kr_VHL_binds_HIF_1_p1)
Rule('VHL_ubiqu_HIF_1_p1',HIF1(p1=('oh', 1)) % VHL(hif_p=1) >> HIF1(p1 = 'ub') + VHL(hif_p=None),
     k_VHL_ubiq_HIF_1_p1)

#HIF degrades after one ubiquination at p1
Rule('HIF_1_p1_binds_proteosome', HIF1(p1='ub') + proteosome(hif_p=None) | HIF1(p1 =('ub', 1)) % VHL(hif_p=1),
    kf_HIF_1_p1_proteo_binding, kr_HIF_1_p1_proteo_binding)
Rule('HIF_1_p1_degrades', HIF1(p1 =('ub', 1)) % VHL(hif_p=1) >> VHL(hif_p=None),
    k_HIF_1_p1_degrades)

#HIF degrades after one ubiquination at p2
Rule('HIF_1_p2_binds_proteosome', HIF1(p2='ub') + proteosome(hif_p=None) | HIF1(p2 =('ub', 1)) % VHL(hif_p=1),
    kf_HIF_1_p1_proteo_binding, kr_HIF_1_p2_proteo_binding)
Rule('HIF_1_p2_degrades', HIF1(p2 =('ub', 1)) % VHL(hif_p=1) >> VHL(hif_p=None),
    k_HIF_1_p1_degrades)

#FIH hydroxylates HIF-1 at a1 site
Rule('HIF_1_binds_FIH', HIF1(a1 = 'u') + FIH(hif_a1=None) | HIF1(a1 =('u', 1)) % FIH(hif_a1=1),
     kf_HIF1_binds_FIH, kr_HIF1_binds_FIH)
Rule('HIF_1_hydroxylated_at_a1', HIF1(a1 =('u', 1)) % FIH(hif_a1=1) >> HIF1(a1 = 'oh') + FIH(hif_a1 = None),
     k_HIF_1_hydroxylated_at_a1)

#p300 interacts with HIF-1
Rule('HIF1_binds_p300', HIF1(a1 = 'u') + p300(hif_a1=None) | HIF1(a1 =('u', 1)) % p300(hif_a1=1),
     kf_HIF1_binds_p300, kr_HIF1_binds_p300)

# Observables

# Observable('degradation', proteosome(hif_p=None))
Observable('p300 binding', p300(hif_a1=None))


# Simulation Commands

tspan = np.linspace(0, 1, 100)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
result = sim.run()

for obs in model.observables:
    plt.plot(tspan, result.observables[obs.name], lw=2, label=obs.name)

plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc=0)

plt.show()


