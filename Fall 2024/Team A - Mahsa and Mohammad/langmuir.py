from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# MODEL #1: Langmuir binding model

Monomer('P', ['b'])  # protein
Monomer('S', ['b']) # surface site

Parameter('P_init', 1000)  # this will probably be in concentration units
Parameter('S_init', 500)

Initial(P(b=None), P_init)
Initial(S(b=None), S_init)

# model dependence of binding on pH as a linear function
Parameter('m', 1)
Parameter('b', 0.1)
Parameter('pH', 7)
Expression('ka', m * pH + b)
Parameter('kd', 1000)
Rule('protein_binds_surface', P(b=None) + S(b=None) | P(b=1) % S(b=1), ka, kd)

Observable('free_protein', P(b=None))
Observable('free_surface_sites', S(b=None))
Observable('bound_surface_sites', P(b=1) % S(b=1))

if __name__ == '__main__':

    # simulations
    PLOT_TIMECOURSES = False
    fraction_sites_bound = []
    ph_vals = np.arange(6, 10.1, 0.2)
    for ph in ph_vals:
        print('pH: %g' % ph)

        tspan = np.linspace(0, 0.002, 101)
        sim = ScipyOdeSimulator(model, tspan, verbose=False)
        output = sim.run(param_values={'pH': ph})

        fraction_sites_bound.append(output.observables['bound_surface_sites'][-1] / S_init.value)

        if PLOT_TIMECOURSES:
            plt.figure(constrained_layout=True)
            plt.title('pH = %g' % ph)
            for obs in model.observables:
                plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
            plt.xlabel('time')
            plt.ylabel('amount')
            plt.legend(loc='best')

    plt.figure(constrained_layout=True)
    plt.plot(ph_vals, fraction_sites_bound, 'o-', ms=8)
    plt.xlabel('pH')
    plt.ylabel('fraction sites bound')

    plt.show()
