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
Monomer('HIF1_nucleus', ['p1','p2', 'a1','gene'],
        {'p1': ['u', 'oh', 'ub'], 'p2': ['u', 'oh', 'ub'], 'a1': ['u', 'oh']})
Monomer('HIF2_nucleus', ['p1','p2', 'a1','gene'],
        {'p1': ['u', 'oh', 'ub'], 'p2': ['u', 'oh', 'ub'], 'a1': ['u', 'oh']})
Monomer('Importer', ['hif'])
Monomer('PHD3_nucleus', ['hif_p'])
Monomer('PHD3_gene', ['prom'])


# Initial Conditions
Initial(HIF1(p1='u', p2='u', a1='u'), Parameter('HIF1_0', 50))
Initial(HIF2(p1='u', p2='u', a1='u'), Parameter('HIF2_0', 50))
Initial(PHD1(hif_p=None), Parameter('PHD1_0', 50))
Initial(PHD2(hif_p=None), Parameter('PHD2_0', 50))
Initial(PHD3(hif_p=None), Parameter('PHD3_0', 50))
Initial(FIH(hif_a1=None), Parameter('FIH_0', 50))
Initial(VHL(hif_p=None), Parameter('VHL_0', 50))
Initial(p300(hif_a1=None), Parameter('p300_0', 50))
Initial(proteosome(hif_p=None), Parameter('proteosome_0', 50))


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

Parameter('kf_HIF1_binds_importer')
Parameter('kr_HIF1_binds_importer')
Parameter('k_HIF1_enters_nucleus')

Parameter('kf_PHD3_binds_HIF1_p2_n')
Parameter('kr_PHD3_binds_HIF1_p2_n')
Parameter('k_PHD3_HIF1_hydroxy_p2_n')

Parameter('kf_PHD3_binds_HIF1_p1_n')
Parameter('kr_PHD3_binds_HIF1_p1_n')
Parameter('k_PHD3_HIF1_hydroxy_p1_n')

Parameter('kf_PHD1_binds_HIF1_p2_n')
Parameter('kr_PHD1_binds_HIF1_p2_n')
Parameter('k_PHD1_HIF1_hydroxy_p2_n')

Parameter('kf_PHD1_binds_HIF1_p1_n')
Parameter('kr_PHD1_binds_HIF1_p1_n')
Parameter('k_PHD1_HIF1_hydroxy_p1_n')

Parameter('kf_HIF1_binds_PHD3_gene')
Parameter('kr_HIF1_binds_PHD3_gene')

Parameter('kf_VHL_binds_HIF_1_p2_n')
Parameter('kr_VHL_binds_HIF_1_p2_n')
Parameter('k_VHL_ubiq_HIF_1_p2_n')

Parameter('kf_VHL_binds_HIF_1_p1_n')
Parameter('kr_VHL_binds_HIF_1_p1_n')
Parameter('k_VHL_ubiq_HIF_1_p1_n')

Parameter('kf_HIF_1_p1_proteo_binding_n')
Parameter('kr_HIF_1_p1_proteo_binding_n')
Parameter('k_HIF_1_p1_degrades_n')

Parameter('kf_HIF_1_p2_proteo_binding_n')
Parameter('kr_HIF_1_p2_proteo_binding_n')
Parameter('k_HIF_1_p2_degrades_n')

# Rules

#PHD3 hydroxylates HIF-1a at Proline 564
Rule('PHD3_binds_HIF1_p2', HIF1(p2='u') + PHD3(hif_p=None) | HIF1(p2=('u', 1)) % PHD3(hif_p=1),
     kf_PHD3_binds_HIF1_p2, kr_PHD3_binds_HIF1_p2 )
Rule('PHD3_hydroxy_HIF1_p2', HIF1(p2=('u', 1)) % PHD3(hif_p=1) >> HIF1(p2='oh') + PHD3(hif_p=None),
     k_PHD3_HIF1_hydroxy_p2)

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
    kf_HIF_1_p2_proteo_binding, kr_HIF_1_p2_proteo_binding)
Rule('HIF_1_p2_degrades', HIF1(p2 =('ub', 1)) % VHL(hif_p=1) >> VHL(hif_p=None),
    k_HIF_1_p2_degrades)

#FIH hydroxylates HIF-1 at a1 site
Rule('HIF_1_binds_FIH', HIF1(a1 = 'u') + FIH(hif_a1=None) | HIF1(a1=('u', 1)) % FIH(hif_a1=1),
     kf_HIF1_binds_FIH, kr_HIF1_binds_FIH)
Rule('HIF_1_hydroxylated_at_a1', HIF1(a1 =('u', 1)) % FIH(hif_a1=1) >> HIF1(a1 = 'oh') + FIH(hif_a1 = None),
     k_HIF_1_hydroxylated_at_a1)

#Representing the HIF-1 protein entering the nucleus
Rule('HIF_1_binds_importer', HIF1(a1 = 'u') + Importer(hif=None) | HIF1(a1=('u', 1)) % Importer(hif=1),
     kf_HIF1_binds_importer, kr_HIF1_binds_importer)
Rule('HIF_1_enters_nucleus',HIF1(a1=('u', 1)) % Importer(hif=1) >> Importer(hif=None) + HIF1_nucleus(a1='u'),
     k_HIF1_enters_nucleus)

#p300 interacts with HIF-1_nucleus
Rule('HIF1_binds_p300', HIF1_nucleus(a1 = 'u') + p300(hif_a1=None) | HIF1_nucleus(a1 =('u', 1)) % p300(hif_a1=1),
     kf_HIF1_binds_p300, kr_HIF1_binds_p300)


#PHD3 hydroxylates HIF-1a at Proline 531
Rule('PHD3_binds_HIF1_p2_in_nucleus', HIF1_nucleus(p2='u') + PHD3_nucleus(hif_p=None) | HIF1_nucleus(p2=('u', 1)) % PHD3_nucleus(hif_p=1),
     kf_PHD3_binds_HIF1_p2_n, kr_PHD3_binds_HIF1_p2_n)
Rule('PHD3_hydroxy_HIF1_p2_in_nucleus', HIF1_nucleus(p2=('u', 1)) % PHD3_nucleus(hif_p=1) >> HIF1_nucleus(p2='oh') + PHD3_nucleus(hif_p=None),
     k_PHD3_HIF1_hydroxy_p2_n)


#PHD1 hydroxylates HIF-1a at Proline 531
Rule('PHD1_binds_HIF1_p2_in_nucleus', HIF1_nucleus(p2='u') + PHD1(hif_p=None) | HIF1_nucleus(p2=('u', 1)) % PHD1(hif_p=1),
     kf_PHD1_binds_HIF1_p2_n, kr_PHD1_binds_HIF1_p2_n)
Rule('PHD1_hydroxy_HIF1_p2_in_nucleus', HIF1(p2=('u', 1)) % PHD1(hif_p=1) >> HIF1(p2='oh') + PHD1(hif_p=None),
     k_PHD1_HIF1_hydroxy_p2_n)

#PHD1 hydroxylated HIF-1a at Proline 405
Rule('PHD1_binds_HIF1_p1_in_nucleus', HIF1_nucleus(p1='u') + PHD1(hif_p=None) | HIF1_nucleus(p1=('u', 1)) % PHD1(hif_p=1),
     kf_PHD1_binds_HIF1_p1_n, kr_PHD1_binds_HIF1_p1_n)
Rule('PHD1_hydroxy_HIF1_p1_in_nucleus', HIF1(p1=('u', 1)) % PHD1(hif_p=1) >> HIF1(p1='oh') + PHD1(hif_p=None),
     k_PHD1_HIF1_hydroxy_p1_n)

#HIF-1 binding to PHD3 gene
Rule('HIF1_binding_to_gene', HIF1_nucleus(a1 =('u', 1), gene = None) % p300(hif_a1=1) + PHD3_gene(prom=None) |
     HIF1_nucleus(a1 =('u', 1), gene=1) % p300(hif_a1=1) % PHD3_gene(prom=1),
     kf_HIF1_binds_PHD3_gene, kr_HIF1_binds_PHD3_gene)

#VHL binds to HIF-1 with 1-OH on p2 and ubiquinates
Rule('VHL_binds_HIF_1_p2_n', HIF1_nucleus(p2='oh') + VHL(hif_p=None) | HIF1_nucleus(p2=('oh', 1)) % VHL(hif_p=1),
     kf_VHL_binds_HIF_1_p2_n, kr_VHL_binds_HIF_1_p2_n)
Rule('VHL_ubiqu_HIF_1_p2_n',HIF1_nucleus(p2=('oh', 1)) % VHL(hif_p=1) >> HIF1_nucleus(p2 = 'ub') + VHL(hif_p=None),
     k_VHL_ubiq_HIF_1_p2_n)

#VHL binds to HIF-1 with 1-OH on p1 and ubiquinates
Rule('VHL_binds_HIF_1_p1_n', HIF1_nucleus(p1='oh') + VHL(hif_p=None) | HIF1_nucleus(p1=('oh', 1)) % VHL(hif_p=1),
     kf_VHL_binds_HIF_1_p1_n, kr_VHL_binds_HIF_1_p1_n)
Rule('VHL_ubiqu_HIF_1_p1_n',HIF1_nucleus(p1=('oh', 1)) % VHL(hif_p=1) >> HIF1_nucleus(p1 = 'ub') + VHL(hif_p=None),
     k_VHL_ubiq_HIF_1_p1_n)

#HIF degrades after one ubiquination at p1
Rule('HIF_1_p1_binds_proteosome_n', HIF1_nucleus(p1='ub') + proteosome(hif_p=None) | HIF1_nucleus(p1 =('ub', 1)) % VHL(hif_p=1),
    kf_HIF_1_p1_proteo_binding_n, kr_HIF_1_p1_proteo_binding_n)
Rule('HIF_1_p1_degrades_n', HIF1_nucleus(p1 =('ub', 1)) % VHL(hif_p=1) >> VHL(hif_p=None),
    k_HIF_1_p1_degrades_n)

#HIF degrades after one ubiquination at p2
Rule('HIF_1_p2_binds_proteosome_n', HIF1_nucleus(p2='ub') + proteosome(hif_p=None) | HIF1_nucleus(p2 =('ub', 1)) % VHL(hif_p=1),
    kf_HIF_1_p2_proteo_binding_n, kr_HIF_1_p2_proteo_binding_n)
Rule('HIF_1_p2_degrades_n', HIF1_nucleus(p2 =('ub', 1)) % VHL(hif_p=1) >> VHL(hif_p=None),
    k_HIF_1_p2_degrades_n)

# Observables
Observable('free_HIF1', HIF1(p1='u', p2='u', a1='u'))
Observable('degradation', proteosome(hif_p=None))
Observable('p300_binding', p300(hif_a1=None))


# Simulation Commands

tspan = np.linspace(0, 30, 301)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
result = sim.run()

for obs in model.observables:
    print('plotting')
    plt.plot(tspan, result.observables[obs.name], lw=2, label=obs.name)

plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc=0)

plt.show()


