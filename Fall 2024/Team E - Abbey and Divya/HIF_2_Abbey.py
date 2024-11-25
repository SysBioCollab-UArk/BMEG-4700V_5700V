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
Parameter('kf_PHD3_binds_HIF2_p2')
Parameter('kr_PHD3_binds_HIF2_p2')
Parameter('k_PHD3_HIF2_hydroxy_p2')

Parameter('kf_PHD3_binds_HIF2_p1')
Parameter('kr_PHD3_binds_HIF2_p1')
Parameter('k_PHD3_HIF2_hydroxy_p1')

Parameter('kf_PHD2_binds_HIF2_p2')
Parameter('kr_PHD2_binds_HIF2_p2')
Parameter('k_PHD2_HIF2_hydroxy_p2')

Parameter('kf_PHD2_binds_HIF2_p1')
Parameter('kr_PHD2_binds_HIF2_p1')
Parameter('k_PHD2_HIF2_hydroxy_p1')

Parameter('kf_VHL_binds_HIF_2_p2')
Parameter('kr_VHL_binds_HIF_2_p2')
Parameter('k_VHL_ubiq_HIF_2_p2')

Parameter('kf_VHL_binds_HIF_2_p1')
Parameter('kr_VHL_binds_HIF_2_p1')
Parameter('k_VHL_ubiq_HIF_2_p1')

Parameter('kf_HIF_2_p1_proteo_binding')
Parameter('kr_HIF_2_p1_proteo_binding')
Parameter('k_HIF_2_p1_degrades')

Parameter('kf_HIF_2_p2_proteo_binding')
Parameter('kr_HIF_2_p2_proteo_binding')
Parameter('k_HIF_2_p2_degrades')

Parameter('kf_HIF2_binds_FIH')
Parameter('kr_HIF2_binds_FIH')
Parameter('k_HIF_2_hydroxylated_at_a1')

Parameter('kf_HIF2_binds_p300')
Parameter('kr_HIF2_binds_p300')

# Rules

#PHD3 hydroxylates HIF-2a at Proline 531
Rule('PHD3_binds_HIF2_p2', HIF2(p2='u') + PHD3(hif_p=None) | HIF2(p2=('u', 1)) % PHD3(hif_p=1),
     kf_PHD3_binds_HIF2_p2, kr_PHD3_binds_HIF2_p2 )
Rule('PHD3_hydroxy_HIF2_p2', HIF2(p2=('u', 1)) % PHD3(hif_p=1) >> HIF2(p2='oh') + PHD3(hif_p=None),
     k_PHD3_HIF2_hydroxy_p2)

#PHD3 hydroxylated HIF-2a at Proline 405
Rule('PHD3_binds_HIF2_p1', HIF2(p1='u') + PHD3(hif_p=None) | HIF2(p1=('u', 1)) % PHD3(hif_p=1),
     kf_PHD3_binds_HIF2_p1, kr_PHD3_binds_HIF2_p1)
Rule('PHD3_hydroxy_HIF2_p1', HIF2(p1=('u', 1)) % PHD3(hif_p=1) >> HIF2(p1='oh') + PHD3(hif_p=None),
     k_PHD3_HIF2_hydroxy_p1)

#PHD2 hydroxylates HIF-2a at Proline 531
Rule('PHD2_binds_HIF2_p2', HIF2(p2='u') + PHD2(hif_p=None) | HIF2(p2=('u', 1)) % PHD2(hif_p=1),
     kf_PHD2_binds_HIF2_p2, kr_PHD2_binds_HIF2_p2 )
Rule('PHD2_hydroxy_HIF2_p2', HIF2(p2=('u', 1)) % PHD2(hif_p=1) >> HIF2(p2='oh') + PHD2(hif_p=None),
     k_PHD2_HIF2_hydroxy_p2)

#PHD2 hydroxylates HIF-2a at Proline 405
Rule('PHD2_binds_HIF2_p1', HIF2(p1='u') + PHD2(hif_p=None) | HIF2(p1=('u', 1)) % PHD2(hif_p=1),
     kf_PHD2_binds_HIF2_p1, kr_PHD2_binds_HIF2_p1)
Rule('PHD2_hydroxy_HIF2_p1', HIF2(p1=('u', 1)) % PHD2(hif_p=1) >> HIF2(p1='oh') + PHD2(hif_p=None),
     k_PHD2_HIF2_hydroxy_p1)

#VHL binds to HIF-2 with 1-OH on p2 and ubiquinates
Rule('VHL_binds_HIF_2_p2', HIF2(p2='oh') + VHL(hif_p=None) | HIF2(p2=('oh', 1)) % VHL(hif_p=1),
     kf_VHL_binds_HIF_2_p2, kr_VHL_binds_HIF_2_p2)
Rule('VHL_ubiqu_HIF_2_p2',HIF2(p2=('oh', 1)) % VHL(hif_p=1) >> HIF2(p2 = 'ub') + VHL(hif_p=None),
     k_VHL_ubiq_HIF_2_p2)

#VHL binds to HIF-2 with 1-OH on p1 and ubiquinates
Rule('VHL_binds_HIF_2_p1', HIF2(p1='oh') + VHL(hif_p=None) | HIF2(p1=('oh', 1)) % VHL(hif_p=1),
     kf_VHL_binds_HIF_2_p1, kr_VHL_binds_HIF_2_p1)
Rule('VHL_ubiqu_HIF_2_p1',HIF2(p1=('oh', 1)) % VHL(hif_p=1) >> HIF2(p1 = 'ub') + VHL(hif_p=None),
     k_VHL_ubiq_HIF_2_p1)

#HIF degrades after one ubiquination at p1
Rule('HIF_2_p1_binds_proteosome', HIF2(p1='ub') + proteosome(hif_p=None) | HIF2(p1 =('ub', 1)) % VHL(hif_p=1),
    kf_HIF_2_p1_proteo_binding, kr_HIF_2_p1_proteo_binding)
Rule('HIF_2_p1_degrades', HIF2(p1 =('ub', 1)) % VHL(hif_p=1) >> VHL(hif_p=None),
    k_HIF_2_p1_degrades)

#HIF degrades after one ubiquination at p2
Rule('HIF_2_p2_binds_proteosome', HIF2(p2='ub') + proteosome(hif_p=None) | HIF2(p2 =('ub', 1)) % VHL(hif_p=1),
    kf_HIF_2_p2_proteo_binding, kr_HIF_2_p2_proteo_binding)
Rule('HIF_2_p2_degrades', HIF2(p2 =('ub', 1)) % VHL(hif_p=1) >> VHL(hif_p=None),
    k_HIF_2_p2_degrades)

#FIH hydroxylates HIF-2 at a1 site
Rule('HIF_2_binds_FIH', HIF2(a1 = 'u') + FIH(hif_a1=None) | HIF2(a1 =('u', 1)) % FIH(hif_a1=1),
     kf_HIF2_binds_FIH, kr_HIF2_binds_FIH)
Rule('HIF_2_hydroxylated_at_a1', HIF2(a1 =('u', 1)) % FIH(hif_a1=1) >> HIF2(a1 = 'oh') + FIH(hif_a1 = None),
     k_HIF_2_hydroxylated_at_a1)

#p300 interacts with HIF-2
Rule('HIF2_binds_p300', HIF2(a1 = 'u') + p300(hif_a1=None) | HIF2(a1 =('u', 1)) % p300(hif_a1=1),
     kf_HIF2_binds_p300, kr_HIF2_binds_p300)


