from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers
Monomer('HIF1', ['p1', 'p2', 'a1', 'loc', 'gene'],
        {'p1': ['u', 'oh', 'ub'], 'p2': ['u', 'oh', 'ub'], 'a1': ['u', 'oh'], 'loc': ['cyt', 'nuc']})
Monomer('HIF2', ['p1','p2', 'a1', 'loc', 'gene'],
        {'p1': ['u', 'oh', 'ub'], 'p2': ['u', 'oh', 'ub'], 'a1': ['u', 'oh'], 'loc': ['cyt', 'nuc']})
Monomer('PHD1', ['hif_p'])
Monomer('PHD2', ['hif_p'])
Monomer('PHD3', ['hif_p'])
Monomer('FIH', ['hif_a1'])
Monomer('VHL', ['hif_p'])
Monomer('p300', ['hif_a1'])
Monomer('proteosome', ['hif_p'])
# Monomer('HIF1_nucleus', ['p1','p2', 'a1','gene'],
#         {'p1': ['u', 'oh', 'ub'], 'p2': ['u', 'oh', 'ub'], 'a1': ['u', 'oh']})
# Monomer('HIF2_nucleus', ['p1','p2', 'a1','gene'],
#         {'p1': ['u', 'oh', 'ub'], 'p2': ['u', 'oh', 'ub'], 'a1': ['u', 'oh']})
Monomer('Importer', ['hif'])
Monomer('PHD3_nucleus', ['hif_p'])
Monomer('PHD3_gene', ['prom'])
Monomer('PHD2_gene', ['prom'])


# Initial Conditions
Initial(HIF1(p1='u', p2='u', a1='u', loc='cyt', gene=None), Parameter('HIF1_0', 50))
Initial(HIF2(p1='u', p2='u', a1='u', loc='cyt', gene=None), Parameter('HIF2_0', 50))
Initial(PHD1(hif_p=None), Parameter('PHD1_0', 50))
Initial(PHD2(hif_p=None), Parameter('PHD2_0', 50))
Initial(PHD3(hif_p=None), Parameter('PHD3_0', 50))
Initial(FIH(hif_a1=None), Parameter('FIH_0', 50))
Initial(VHL(hif_p=None), Parameter('VHL_0', 50))
Initial(p300(hif_a1=None), Parameter('p300_0', 50))
Initial(proteosome(hif_p=None), Parameter('proteosome_0', 50))
Initial(Importer(hif=None), Parameter('importer_0', 50))
Initial(PHD3_gene(prom=None), Parameter('PHD3_gene_0', 50))
Initial(PHD2_gene(prom=None), Parameter('PHD2_gene_0', 50))

# Parameters for HIF1 pathway
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

Parameter('kf_PHD1_binds_HIF1_p2_n')
Parameter('kr_PHD1_binds_HIF1_p2_n')
Parameter('k_PHD1_HIF1_hydroxy_p2_n')

Parameter('kf_PHD1_binds_HIF1_p1_n')
Parameter('kr_PHD1_binds_HIF1_p1_n')
Parameter('k_PHD1_HIF1_hydroxy_p1_n')

Parameter('kf_HIF1_binds_PHD3_gene')
Parameter('kr_HIF1_binds_PHD3_gene')

Parameter('kf_HIF1_binds_PHD2_gene')
Parameter('kr_HIF1_binds_PHD2_gene')

Parameter('kf_hif1_PHD3_created')
Parameter('kr_hif1_PHD3_created')

Parameter('kf_hif1_PHD2_created')
Parameter('kr_hif1_PHD2_created')

# Parameters for HIF2 Pathway
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

Parameter('kf_HIF2_binds_importer')
Parameter('kr_HIF2_binds_importer')
Parameter('k_HIF2_enters_nucleus')

Parameter('kf_PHD1_binds_HIF2_p2_n')
Parameter('kr_PHD1_binds_HIF2_p2_n')
Parameter('k_PHD1_HIF2_hydroxy_p2_n')

Parameter('kf_PHD1_binds_HIF2_p1_n')
Parameter('kr_PHD1_binds_HIF2_p1_n')
Parameter('k_PHD1_HIF2_hydroxy_p1_n')

Parameter('kf_HIF2_binds_PHD3_gene')
Parameter('kr_HIF2_binds_PHD3_gene')

Parameter('kf_hif2_PHD3_created')
Parameter('kr_hif2_PHD3_created')

# Rules for HIF1

#PHD3 hydroxylates HIF-1a at Proline 564
Rule('PHD3_binds_HIF1_p2', HIF1(p2='u') + PHD3(hif_p=None) | HIF1(p2=('u', 1)) % PHD3(hif_p=1),
     kf_PHD3_binds_HIF1_p2, kr_PHD3_binds_HIF1_p2 )
Rule('PHD3_hydroxy_HIF1_p2', HIF1(p2=('u', 1)) % PHD3(hif_p=1) >> HIF1(p2='oh') + PHD3(hif_p=None),
     k_PHD3_HIF1_hydroxy_p2)

#PHD2 hydroxylates HIF-1a at Proline 564
Rule('PHD2_binds_HIF1_p2', HIF1(p2='u', loc='cyt') + PHD2(hif_p=None) | HIF1(p2=('u', 1), loc='cyt') % PHD2(hif_p=1),
     kf_PHD2_binds_HIF1_p2, kr_PHD2_binds_HIF1_p2 )
Rule('PHD2_hydroxy_HIF1_p2', HIF1(p2=('u', 1), loc='cyt') % PHD2(hif_p=1) >> HIF1(p2='oh', loc='cyt') + PHD2(hif_p=None),
     k_PHD2_HIF1_hydroxy_p2)

#PHD2 hydroxylates HIF-1a at Proline 402
Rule('PHD2_binds_HIF1_p1', HIF1(p1='u', loc='cyt') + PHD2(hif_p=None) | HIF1(p1=('u', 1), loc='cyt') % PHD2(hif_p=1),
     kf_PHD2_binds_HIF1_p1, kr_PHD2_binds_HIF1_p1)
Rule('PHD2_hydroxy_HIF1_p1', HIF1(p1=('u', 1), loc='cyt') % PHD2(hif_p=1) >> HIF1(p1='oh', loc='cyt') + PHD2(hif_p=None),
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
Rule('HIF_1_binds_FIH', HIF1(a1 = 'u', loc='cyt') + FIH(hif_a1=None) | HIF1(a1=('u', 1), loc='cyt') % FIH(hif_a1=1),
     kf_HIF1_binds_FIH, kr_HIF1_binds_FIH)
Rule('HIF_1_hydroxylated_at_a1', HIF1(a1 =('u', 1), loc='cyt') % FIH(hif_a1=1) >> HIF1(a1 = 'oh', loc='cyt') + FIH(hif_a1 = None),
     k_HIF_1_hydroxylated_at_a1)

#Representing the HIF-1 protein entering the nucleus
Rule('HIF_1_binds_importer', HIF1(a1 = 'u', loc='cyt') + Importer(hif=None) | HIF1(a1=('u', 1), loc='cyt') % Importer(hif=1),
     kf_HIF1_binds_importer, kr_HIF1_binds_importer)
Rule('HIF_1_enters_nucleus',HIF1(a1=('u', 1), loc='cyt') % Importer(hif=1) >> Importer(hif=None) + HIF1(a1='u', loc='nuc'),
     k_HIF1_enters_nucleus)

#p300 interacts with HIF-1_nucleus
Rule('HIF1_binds_p300', HIF1(a1 = 'u', loc='nuc') + p300(hif_a1=None) | HIF1(a1 =('u', 1), loc='nuc') % p300(hif_a1=1),
     kf_HIF1_binds_p300, kr_HIF1_binds_p300)

#PHD1 hydroxylates HIF-1a at Proline 564
Rule('PHD1_binds_HIF1_p2_in_nucleus', HIF1(p2='u', loc='nuc') + PHD1(hif_p=None) | HIF1(p2=('u', 1), loc='nuc') % PHD1(hif_p=1),
     kf_PHD1_binds_HIF1_p2_n, kr_PHD1_binds_HIF1_p2_n)
Rule('PHD1_hydroxy_HIF1_p2_in_nucleus', HIF1(p2=('u', 1), loc='nuc') % PHD1(hif_p=1) >> HIF1(p2='oh', loc='nuc') + PHD1(hif_p=None),
     k_PHD1_HIF1_hydroxy_p2_n)

#PHD1 hydroxylated HIF-1a at Proline 402
Rule('PHD1_binds_HIF1_p1_in_nucleus', HIF1(p1='u', loc='nuc') + PHD1(hif_p=None) | HIF1(p1=('u', 1), loc='nuc') % PHD1(hif_p=1),
     kf_PHD1_binds_HIF1_p1_n, kr_PHD1_binds_HIF1_p1_n)
Rule('PHD1_hydroxy_HIF1_p1_in_nucleus', HIF1(p1=('u', 1), loc='nuc') % PHD1(hif_p=1) >> HIF1(p1='oh', loc='nuc') + PHD1(hif_p=None),
     k_PHD1_HIF1_hydroxy_p1_n)

#HIF-1 binding to PHD3 gene
Rule('HIF1_binding_to_gene_PHD3', HIF1(a1 =('u', 1), gene = None, loc='nuc') % p300(hif_a1=1) + PHD3_gene(prom=None) |
     HIF1(a1 =('u', 1), gene=2, loc='nuc') % p300(hif_a1=1) % PHD3_gene(prom=2),
     kf_HIF1_binds_PHD3_gene, kr_HIF1_binds_PHD3_gene)

#HIF-1 binding to PHD2 gene
Rule('HIF1_binding_to_gene_PHD2', HIF1(a1 =('u', 1), gene = None, loc='nuc') % p300(hif_a1=1) + PHD2_gene(prom=None) |
     HIF1(a1 =('u', 1), gene=2, loc='nuc') % p300(hif_a1=1) % PHD2_gene(prom=2),
     kf_HIF1_binds_PHD2_gene, kr_HIF1_binds_PHD2_gene)

# HIF1 - PHD3 gene making a protein
Rule('HIF1_making_PHD3', HIF1(a1=('u', 1), gene=2, loc='nuc') % p300(hif_a1=1) % PHD3_gene(prom=2) |
     HIF1(a1=('u', 1), gene=2, loc='nuc') % p300(hif_a1=1) % PHD3_gene(prom=2) + PHD3(hif_p=None),
     kf_hif1_PHD3_created, kr_hif1_PHD3_created)

# HIF1 - PHD2 gene making a protein
Rule('HIF1_making_PHD2', HIF1(a1=('u', 1), gene=2, loc='nuc') % p300(hif_a1=1) % PHD2_gene(prom=2) |
     HIF1(a1=('u', 1), gene=2, loc='nuc') % p300(hif_a1=1) % PHD2_gene(prom=2) + PHD2(hif_p=None),
     kf_hif1_PHD2_created, kr_hif1_PHD2_created)



# Rules for HIF2
#PHD3 hydroxylates HIF-2a at Proline 531
Rule('PHD3_binds_HIF2_p2', HIF2(p2='u') + PHD3(hif_p=None) | HIF2(p2=('u', 1)) % PHD3(hif_p=1),
     kf_PHD3_binds_HIF2_p2, kr_PHD3_binds_HIF2_p2 )
Rule('PHD3_hydroxy_HIF2_p2', HIF2(p2=('u', 1)) % PHD3(hif_p=1) >> HIF2(p2='oh') + PHD3(hif_p=None),
     k_PHD3_HIF2_hydroxy_p2)

#PHD2 hydroxylates HIF-2a at Proline 531
Rule('PHD2_binds_HIF2_p2', HIF2(p2='u', loc='cyt') + PHD2(hif_p=None) | HIF2(p2=('u', 1), loc='cyt') % PHD2(hif_p=1),
     kf_PHD2_binds_HIF2_p2, kr_PHD2_binds_HIF2_p2 )
Rule('PHD2_hydroxy_HIF2_p2', HIF2(p2=('u', 1), loc='cyt') % PHD2(hif_p=1) >> HIF2(p2='oh', loc='cyt') + PHD2(hif_p=None),
     k_PHD2_HIF2_hydroxy_p2)

#PHD2 hydroxylates HIF-2a at Proline 405
Rule('PHD2_binds_HIF2_p1', HIF2(p1='u', loc='cyt') + PHD2(hif_p=None) | HIF2(p1=('u', 1), loc='cyt') % PHD2(hif_p=1),
     kf_PHD2_binds_HIF2_p1, kr_PHD2_binds_HIF2_p1)
Rule('PHD2_hydroxy_HIF2_p1', HIF2(p1=('u', 1), loc='cyt') % PHD2(hif_p=1) >> HIF2(p1='oh', loc='cyt') + PHD2(hif_p=None),
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
Rule('HIF_2_binds_FIH', HIF2(a1 = 'u', loc='cyt') + FIH(hif_a1=None) | HIF2(a1=('u', 1), loc='cyt') % FIH(hif_a1=1),
     kf_HIF2_binds_FIH, kr_HIF2_binds_FIH)
Rule('HIF_2_hydroxylated_at_a1', HIF2(a1 =('u', 1), loc='cyt') % FIH(hif_a1=1) >> HIF2(a1 = 'oh', loc='cyt') + FIH(hif_a1 = None),
     k_HIF_2_hydroxylated_at_a1)

#Representing the HIF-2 protein entering the nucleus
Rule('HIF_2_binds_importer', HIF2(a1 = 'u', loc='cyt') + Importer(hif=None) | HIF2(a1=('u', 1), loc='cyt') % Importer(hif=1),
     kf_HIF2_binds_importer, kr_HIF2_binds_importer)
Rule('HIF_2_enters_nucleus',HIF2(a1=('u', 1), loc='cyt') % Importer(hif=1) >> Importer(hif=None) + HIF2(a1='u', loc='nuc'),
     k_HIF2_enters_nucleus)

#p300 interacts with HIF-2_nucleus
Rule('HIF2_binds_p300', HIF2(a1 = 'u', loc='nuc') + p300(hif_a1=None) | HIF2(a1 =('u', 1), loc='nuc') % p300(hif_a1=1),
     kf_HIF2_binds_p300, kr_HIF2_binds_p300)

#PHD1 hydroxylates HIF-2a at Proline 531
Rule('PHD1_binds_HIF2_p2_in_nucleus', HIF2(p2='u', loc='nuc') + PHD1(hif_p=None) | HIF2(p2=('u', 1), loc='nuc') % PHD1(hif_p=1),
     kf_PHD1_binds_HIF2_p2_n, kr_PHD1_binds_HIF2_p2_n)
Rule('PHD1_hydroxy_HIF2_p2_in_nucleus', HIF2(p2=('u', 1), loc='nuc') % PHD1(hif_p=1) >> HIF2(p2='oh', loc='nuc') + PHD1(hif_p=None),
     k_PHD1_HIF2_hydroxy_p2_n)

#PHD1 hydroxylated HIF-2a at Proline 405
Rule('PHD1_binds_HIF2_p1_in_nucleus', HIF2(p1='u', loc='nuc') + PHD1(hif_p=None) | HIF2(p1=('u', 1), loc='nuc') % PHD1(hif_p=1),
     kf_PHD1_binds_HIF2_p1_n, kr_PHD1_binds_HIF2_p1_n)
Rule('PHD1_hydroxy_HIF2_p1_in_nucleus', HIF2(p1=('u', 1), loc='nuc') % PHD1(hif_p=1) >> HIF2(p1='oh', loc='nuc') + PHD1(hif_p=None),
     k_PHD1_HIF2_hydroxy_p1_n)

#HIF-2 binding to PHD3 gene
Rule('HIF2_binding_to_gene', HIF2(a1=('u', 1), gene=None, loc='nuc') % p300(hif_a1=1) + PHD3_gene(prom=None) |
     HIF2(a1=('u', 1), gene=2, loc='nuc') % p300(hif_a1=1) % PHD3_gene(prom=2),
     kf_HIF2_binds_PHD3_gene, kr_HIF2_binds_PHD3_gene)

# HIF2 - PHD3 gene making a protein
Rule('HIF2_making_PHD3', HIF2(a1=('u', 1), gene=2, loc='nuc') % p300(hif_a1=1) % PHD3_gene(prom=2) |
     HIF2(a1=('u', 1), gene=2, loc='nuc') % p300(hif_a1=1) % PHD3_gene(prom=2) + PHD3(hif_p=None),
     kf_hif2_PHD3_created, kr_hif2_PHD3_created)


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


