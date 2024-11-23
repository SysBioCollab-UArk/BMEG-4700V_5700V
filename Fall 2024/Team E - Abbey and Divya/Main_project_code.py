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

Parameter('kf_PHD2_binds_HIF1')
Parameter('kr_PHD2_binds_HIF1')
Parameter('k_PHD2_HIF1_hydroxy')

Rule('PHD2_binds_HIF1', HIF1(p1='u') + PHD2(hif_p=None) | HIF1(p1=('u', 1)) % PHD2(hif_p=1),
     kf_PHD2_binds_HIF1, kr_PHD2_binds_HIF1)
Rule('PHD2_HIF1_hydroxy', HIF1(p1=('u', 1)) % PHD2(hif_p=1) >> HIF1(p1='oh') + PHD2(hif_p=None),
     k_PHD2_HIF1_hydroxy)

# Initial Conditions
#
# Parameter('A_0', 100)
# Parameter('B_0', 80)
# Initial(A(b=None), A_0)
# Initial(B(a=None, state='u'), B_0)
