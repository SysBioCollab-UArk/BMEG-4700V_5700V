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

#Initializing Variables
Monomer("ECM", ["itgb1"])
Monomer("ITGB1", ["ecm"])
Monomer("FAK", ["pi3k","state"], {"state": ["u","p"]})
Monomer("PI3K",["fak","state"], {"state": ["i","a"]})
Monomer("PIP", ["state"], {"state": ["_2","_3"]})
Monomer("AKT", ["state"], {"state": ["u","p"]})
Monomer("METABOLIC", ["state"], {"state": ["i","a"]})

#Initials
Initial(ECM(itgb1=None), Parameter('ECM_0', 50))
Initial(ITGB1(ecm=None), Parameter('ITGB1_0', 150))
Initial(FAK(pi3k=None, state='u'), Parameter('FAK_0', 30))
Initial(PI3K(fak=None, state='i'), Parameter('PI3K_0', 40))
Initial(PIP(state='_2'), Parameter('PIP_0', 50))
Initial(AKT(state='u'), Parameter('AKT_0', 70))
Initial(METABOLIC(state='i'), Parameter('met_0', 500))


#Integrin-b1 binds to the ECM (reversible)
Parameter('kf_itgb1_ecm', 0.1)
#Parameter('kr_itgb1_ecm', 10)
Rule("ITGB1_binds_ECM", ITGB1(ecm=None) + ECM(itgb1=None) >> ITGB1(ecm=1) % ECM(itgb1=1),
     kf_itgb1_ecm)

#Integrin-b1 phosphorylates FAK.
Parameter('k_FAK_p', 0.02)
Rule("ITGB1_ECM_phos_FAK", ITGB1(ecm=1) % ECM(itgb1=1) + FAK(pi3k=None, state="u") >>
     ITGB1(ecm=None) + ECM(itgb1=None) + FAK(pi3k=None, state="p"), k_FAK_p)

#Phosphorylated FAK then binds to PI3K
Parameter('kf_FAKp_PI3Ki', 0.1)
Parameter('kr_FAKP_PI3Ki', 0.2)
Rule("FAKp_binds_PI3Ki", FAK(pi3k=None, state="p") + PI3K(fak=None, state="i") |
     FAK(pi3k=1, state="p") % PI3K(fak=1, state="i"), kf_FAKp_PI3Ki, kr_FAKP_PI3Ki)

#The inactive PI3K becomes active
Parameter('k_PI3Ki_a', 0.4)
Rule("PI3Ki_activates", FAK(pi3k=1, state="p") % PI3K(fak=1, state="i") >>
     FAK(pi3k=None, state="u") + PI3K(fak=None, state="a"), k_PI3Ki_a)

#The activated PI3K produces PIP3
Parameter('kPI3Ka_PIP3', 0.4)
Rule("PI3Ka_produces_PIP3", PI3K(fak=None, state="a") + PIP(state = "_2") >> PI3K(fak=None, state="i")
     + PIP(state = "_3"), kPI3Ka_PIP3)

#PI3K inactivates
#Parameter('k_PI3Ka_i', 1)
#Rule("PI3Ka_inactivates", PI3K(fak=None, state="a") >> PI3K(fak=None, state="i"), k_PI3Ka_i)

#PIP3 then phosphorylates AKT
Parameter('k_PIP3_AKTp', 1)
Rule("PIP3_phos_AKTu", PIP(state = "_3") + AKT(state = "u") >>  PIP(state = "_2") + AKT(state = "p"), k_PIP3_AKTp)

#AKT then unphosphorylates by activating a separate metabolic chain.
Parameter('k_AKTp_u', 0.1)
Rule("AKTp_unphos", AKT(state = "p") + METABOLIC(state = "i") >> AKT(state = "u")+ METABOLIC(state = "a"), k_AKTp_u)

Observable('ECM_bindings', ITGB1(ecm=1) % ECM(itgb1=1))
Observable('Metabolic_Activations', METABOLIC(state = 'a'))
#Observable('S_free', X(b=None, state='S'))
#Observable('ES_complex', E(b=ANY))
#Observable('Product', X(state='P'))
#Observable('Better_be_zero', X(b=ANY, state='P'))

if __name__ == '__main__':

    # Simulation commands
    tspan = np.linspace(0, 40, 301)
    sim = ScipyOdeSimulator(model, tspan, verbose=True)
    output = sim.run()

    for obs in model.observables:
        plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
    plt.xlabel('time')
    plt.ylabel('# of molecules')
    plt.legend(loc=0)
    plt.tight_layout()

    plt.show()
