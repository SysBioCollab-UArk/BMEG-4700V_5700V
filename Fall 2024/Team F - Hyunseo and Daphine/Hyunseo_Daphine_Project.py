from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

print(model)
# Roles:
# •	Integrin-β1: Prominent in cell-cell adhesion and ECM interactions.
# •	FAK (Focal Adhesion Kinase):
# •	PI3K (Phosphoinositide 3-kinase):
# •	AKT (Protein Kinase B):
# Detailed Mechanism:
# 1.	ITGb1 attaches to ECM structural protein collagen I.
# 2.	The attachment calls for a conformational change in ITGb1 that activates FAK by autophosphorylation.
# 3.	The activated FAK then can bind to one of the subunits of PI3K activating the enzyme to induce production of second messenger PIP3 (Phosphatidylinositol (3,4,5)-triphosphate).
# 4.	PIP3 then activates AKT and then phosphorylated AKT (pAKT) then can phosphorylate various of other messengers to either promote or inhibit different cellular activity.

Monomer("ECM", ["itgb1"])
Monomer("ITGB1", ["ecm"])
Monomer("FAK", ["pi3k","state"], {"state": ["u","p"]})
Monomer("PI3K",["fak","state"], {"state": ["i","a"]})
Monomer("PIP", ["state"], {"state": ["_2","_3"]})
Monomer("AKT", ["state"], {"state": ["u","p"]})

#
Rule("ITGB1_binds_ECM", ITGB1(ecm=None) + ECM(itgb1=None) | ITGB1(ecm=1) % ECM(itgb1=1),
     kf_itgb1_ecm, kr_itgb1_ecm)

#
Rule("ITGB1_ECM_phos_FAK", ITGB1(ecm=1) % ECM(itgb1=1) + FAK(pi3k=None, state="u") >>
     ITGB1(ecm=1) % ECM(itgb1=1) + FAK(pi3k=None, state="p"), k_FAK_p)

#
Rule("FAKp_binds_PI3Ki", FAK(pi3k=None, state="p") + PI3K(fak=None, state="i") |
     FAK(pi3k=1, state="p") % PI3K(fak=1, state="i"), kf, kr)

Rule("PI3Ki_activates", FAK(pi3k=1, state="p") % PI3K(fak=1, state="i") >>
     FAK(pi3k=None, state="u") + PI3K(fak=None, state="a"), k)

Rule("PI3Ka_produces_PIP3", PI3K(fak=None, state="a") + PIP(state = "_2") >> PI3K(fak=None, state="a")
     + PIP(state = "_3"), k1)

Rule("PI3Ka_inactivates", PI3K(fak=None, state="a") >> PI3K(fak=None, state="i"), k2)

Rule("PIP3_phos_AKTu", PIP(state = "_3") + AKT(state = "u") >>  PIP(state = "_2") + AKT(state = "p"), k3)

Rule("AKTp_unphos", AKT(state = "p") >> AKT(state = "u"), k4)