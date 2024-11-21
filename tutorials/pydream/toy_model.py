from pysb import *

# Empty model
Model()

# Define Monomers
# Enzyme
Monomer('E', ['binding1'])
# Substrate 1
Monomer('S1', ['binding', 'state'], {'state': ['sub', 'pro']})
# Substrate 2
Monomer('S2', ['binding', 'state'], {'state': ['sub', 'pro']})

# Initial Conditions
Parameter('E_0', 1000)
Parameter('S1_0', 500)
Parameter('S2_0', 600)
Initial(E(binding1=None), E_0)
Initial(S1(binding=None, state='sub'), S1_0)
Initial(S2(binding=None, state='sub'), S2_0)

# Observables
Observable('enzyme_total', E())
Observable('enzyme_free', E(binding1=None))
Observable('substrate_1', S1(binding=None, state='sub'))
Observable('substrate_2', S2(binding=None, state='sub'))
Observable('complex_1', E(binding1=1) % S1(binding=1, state='sub'))
Observable('complex_2', E(binding1=1) % S2(binding=1, state='sub'))
Observable('product_1', S1(binding=None, state='pro'))
Observable('product_2', S2(binding=None, state='pro'))
Observable('product_total', S1(state='pro') + S2(state='pro'))

# Define the binding rules
Parameter('k_1', 0.002)
Parameter('k_2', 0.001)
Rule('binding_1',
     E(binding1=None) + S1(state='sub', binding=None) | E(binding1=1) % S1(state='sub', binding=1), k_1, k_2)
Parameter('k_4', 0.004)
#PYDREAM_IT prior k_5 uniform 2
Parameter('k_5', 0.001)
Rule('binding_2',
     E(binding1=None) + S2(state='sub', binding=None) | E(binding1=1) % S2(state='sub', binding=1), k_4, k_5)

# Catalyze
Parameter('k_3', 0.1)
#PYDREAM_IT no-sample k_6
Parameter('k_6', 0.1)
Rule('catalyze_1', E(binding1=1) % S1(state='sub', binding=1) >> E(binding1=None) + S1(state='pro', binding=None), k_3)
Rule('catalyze_2', E(binding1=1) % S2(state='sub', binding=1) >> E(binding1=None) + S2(state='pro', binding=None), k_6)

if __name__ == "__main__":
    from pysb.simulator import ScipyOdeSimulator
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import norm
    import csv

    tspan = np.linspace(0, 30, num=301)  # 601
    solver = ScipyOdeSimulator(model, tspan=tspan, verbose=True)
    output = solver.run()

    plt.figure(figsize=(6.4*1.2, 4.8))
    color = {}
    for obs in model.observables:
        p = plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
        color[obs.name] = p[0].get_color()

    # Create synthetic data
    t_sample_idx = [i for i in range(0, len(tspan), 20)]
    synth_obs = ['complex_2', 'product_total']
    synth_data = {}
    for obs in synth_obs:
        # 'random_state' is the random seed. Set to None if you want to create a different set of data.
        synth_data[obs] = norm.rvs(output.observables[obs][t_sample_idx], 0.1 * output.observables[obs][t_sample_idx],
                                   random_state=10)
        plt.plot(tspan[t_sample_idx], synth_data[obs], 'o', ms=6, color=color[obs])
    plt.xlabel('time')
    plt.ylabel('concentration')
    plt.legend(loc=0, bbox_to_anchor=[1, 1])
    plt.tight_layout()

    # Output data to csv files, if desired
    output_to_file = True  # Set to True if you want to save the synthetic data for calibration
    if output_to_file:
        with open('toy_model_synth_data.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(['observable', 'time', 'time_units', 'average', 'stderr', 'amount_units', 'expt_id'])
            for obs in synth_obs:
                for i,t in enumerate(t_sample_idx):
                    writer.writerow([obs, tspan[t], 'arbitrary', synth_data[obs][i],
                                       max(5, 0.1 * output.observables[obs][t]), 'arbitrary', 'tutorial'])

    plt.show()
