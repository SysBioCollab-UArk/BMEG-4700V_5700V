"""A simple three-species chemical kinetics system known as "Robertson's
example", as presented in:

H. H. Robertson, The solution of a set of reaction rate equations, in Numerical
Analysis: An Introduction, J. Walsh, ed., Academic Press, 1966, pp. 178-182.
"""

# This is a simple system often used to study stiffness in systems of
# differential equations.  It doesn't leverage the power of rules-based modeling
# or pysb, but it's a useful small model for purposes of experimentation.
#
# A brief report addressing issues of stiffness encountered in numerical
# integration of Robertson's example can be found here:
# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.53.8603&rep=rep1&type=pdf
#
# The chemical model is as follows:
#
#      Reaction        Rate
#   ------------------------
#       A -> B         0.04
#      2B -> B + C     3.0e7
#   B + C -> A + C     1.0e4
#
# The resultant system of differential equations is:
#
# y1' = -0.04 * y1 + 1.0e4 * y2 * y3
# y2' =  0.04 * y1 - 1.0e4 * y2 * y3 - 3.0e7 * y2^2
# y3' =                                3.0e7 * y2^2

from pysb import *

Model()

Monomer('A')
Monomer('B')
Monomer('C')

#       A -> B         0.04
Rule('A_to_B', A() >> B(), Parameter('k1', 0.04))

#      2B -> B + C     3.0e7
Rule('BB_to_BC', B() + B() >> B() + C(), Parameter('k2', 3.0e7))

#   B + C -> A + C     1.0e4
Rule('BC_to_AC', B() + C() >> A() + C(), Parameter('k3', 1.0e4))

# The system is known to be stiff for initial values A=1, B=0, C=0
Initial(A(), Parameter('A_0', 1.0))
Initial(B(), Parameter('B_0', 0.0))
Initial(C(), Parameter('C_0', 0.0))

# Observe total amount of each monomer
Observable('A_total', A())
Observable('B_total', B())
Observable('C_total', C())

if __name__ == '__main__':
    from pysb.simulator import ScipyOdeSimulator
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import norm
    import csv

    tspan = np.linspace(0, 40, 401)
    sim = ScipyOdeSimulator(model, tspan, verbose=True)
    output = sim.run()

    color = {}
    for obs in model.observables:
        p = plt.plot(tspan, output.observables[obs.name], label=obs.name)
        color[obs.name] = p[0].get_color()
    plt.xlabel('time')
    plt.ylabel('concentration')
    plt.legend(loc='best')
    plt.tight_layout()

    # Create synthetic data
    t_sample_idx = [i for i in range(0, len(tspan), 40)]
    synth_obs = ['C_total']
    synth_data = {}
    for obs in synth_obs:
        # 'random_state' is the random seed. Set to None if you want to create a different set of data.
        synth_data[obs] = norm.rvs(output.observables[obs][t_sample_idx], 0.35 * output.observables[obs][t_sample_idx],
                                   random_state=10)
        plt.plot(tspan[t_sample_idx], synth_data[obs], 'o', ms=6, color=color[obs])

    # Output synthetic data to csv files, if desired
    output_to_file = True  # Set to True if you want to save the synthetic data for calibration
    if output_to_file:
        with open('robertson_synth_data.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(['observable', 'time', 'time_units', 'average', 'stderr', 'amount_units', 'expt_id'])
            for obs in synth_obs:
                for i, t in enumerate(t_sample_idx):
                    writer.writerow([obs, tspan[t], 'arbitrary', synth_data[obs][i],
                                     max(1e-2, 0.35 * output.observables[obs][t]), 'arbitrary', 'tutorial'])

    plt.show()
