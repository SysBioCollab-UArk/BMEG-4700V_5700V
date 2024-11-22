from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
Model()
# Monomers

Monomer('HIF1', ['p1', 'p2', 'a1'],
        {'p1':['u', 'oh', 'ub'], 'p2':['u', 'oh', 'ub'], 'a1':['u', 'oh']})
Monomer('HIF2', ['p1','p2', 'a1'],
        {'p1':['u', 'oh', 'ub'], 'p2':['u', 'oh', 'ub'], 'a1':['u', 'oh']})
Monomer('PHD1', ['hif_p'])
Monomer('PHD2', ['hif_p'])
Monomer('PHD3', ['hif_p'])
Monomer('FIH', ['hif_a1'])
Monomer('VHL', ['hif_p'])
Monomer('p300', ['hif_a1'])

Parameter('kr_AB')
Parameter('kf_AB')

Rule('PHD2_binds_HIF1', HIF1(p1=('u', None)) + PHD2(hif_p=None) | HIF1(p1=('u', 1)) % PHD(hif_p=1), kf_AB, kr_AB)
Rule('PHD2_HIF1_hydroxy', HIF1(p1=(1, 'u')) % PHD(hif_p=1) >> HIF1(p1=(None, 'oh')) + PHD(hif_p=None), kf_AB)

# Initial Conditions
#
# Parameter('A_0', 100)
# Parameter('B_0', 80)
# Initial(A(b=None), A_0)
# Initial(B(a=None, state='u'), B_0)
