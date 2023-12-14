import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator
from redox_ratio import model
from sklearn.metrics import mean_squared_error

# Load experimental data from different sheets
data_nadhf = pd.read_excel('/mnt/data/elife-73808-fig2-data1-v2.xlsx', sheet_name='nadhf')
data_nadhb = pd.read_excel('/mnt/data/elife-73808-fig2-data1-v2.xlsx', sheet_name='nadhb')
data_bound_ratio = pd.read_excel('/mnt/data/elife-73808-fig2-data1-v2.xlsx', sheet_name='bound ratio')

# Assuming result from optimization (replace with actual optimized parameters)
fitted_params = [100, 50, 0, 0, 0, 0.01, 1, 0.01, 1, 0.01, 1]  # Example values

# Define time span for the simulation
tspan = np.linspace(0, len(data_nadhf) - 1, len(data_nadhf))

# Function to simulate the model
def simulate_model(params, model, tspan):
    for i, param in enumerate(model.parameters):
        if i < len(params):
            param.value = params[i]
    sim = ScipyOdeSimulator(model, tspan)
    return sim.run()

# Simulate the model with optimized parameters
optimized_simulation = simulate_model(fitted_params, model, tspan)
nadh_ratio = optimized_simulation.observables['NADH_free'] / optimized_simulation.observables['NADH_bound']

# Plotting results for each inhibitor condition
conditions = ['oxamate', 'rot', 'oligo', 'fccp']
for condition in conditions:
    plt.figure(figsize=(15, 5))

    # Plot for NADH free
    plt.subplot(1, 3, 1)
    plt.plot(tspan, optimized_simulation.observables['NADH_free'], label='Model NADH_free')
    plt.scatter(np.arange(len(data_nadhf)), data_nadhf[condition], color='red', label=f'Exp NADH_free ({condition})')
    plt.xlabel("Time")
    plt.ylabel("NADH Free")
    plt.legend()

    # Plot for NADH bound
    plt.subplot(1, 3, 2)
    plt.plot(tspan, optimized_simulation.observables['NADH_bound'], label='Model NADH_bound')
    plt.scatter(np.arange(len(data_nadhb)), data_nadhb[condition], color='red', label=f'Exp NADH_bound ({condition})')
    plt.xlabel("Time")
    plt.ylabel("NADH Bound")
    plt.legend()

    # Plot for NADH ratio
    plt.subplot(1, 3, 3)
    plt.plot(tspan, nadh_ratio, label='Model NADH_ratio')
    plt.scatter(np.arange(len(data_bound_ratio)), data_bound_ratio[condition], color='red', label=f'Exp NADH_ratio ({condition})')
    plt.xlabel("Time")
    plt.ylabel("NADH Ratio")
    plt.legend()

    plt.suptitle(f"NADH Dynamics under {condition} Condition")
    plt.tight_layout()
    plt.show()
