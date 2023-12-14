import pandas as pd
import numpy as np
from pysb.simulator import ScipyOdeSimulator
from redox_ratio import model
import json
from scipy.optimize import minimize
# Load experimental data from different sheets
data_nadhf = pd.read_excel('/mnt/data/elife-73808-fig2-data1-v2.xlsx', sheet_name='nadhf')
data_nadhb = pd.read_excel('/mnt/data/elife-73808-fig2-data1-v2.xlsx', sheet_name='nadhb')
data_bound_ratio = pd.read_excel('/mnt/data/elife-73808-fig2-data1-v2.xlsx', sheet_name='bound ratio')

# Define the time span for the simulation
tspan = np.linspace(0, len(data_nadhf) - 1, len(data_nadhf))

# Define the cost function
def cost_function(params, model, tspan, data_nadhf, data_nadhb, data_bound_ratio):
    total_cost = 0

    # Update model parameters
    for i, param in enumerate(model.parameters):
        if i < len(params):
            param.value = params[i]

    # Simulate the model
    sim = ScipyOdeSimulator(model, tspan)
    simulation_results = sim.run()

    # Iterating through each inhibitor condition
    for condition in ['oxamate', 'rot', 'oligo', 'fccp']:
        simulated_nadh_free = simulation_results.observables['NADH_free']
        simulated_nadh_bound = simulation_results.observables['NADH_bound']
        simulated_nadh_ratio = simulation_results.observables['NADH_ratio']

        experimental_nadh_free = data_nadhf[condition].values
        experimental_nadh_bound = data_nadhb[condition].values
        experimental_nadh_ratio = data_bound_ratio[condition].values

        # Calculate cost for each condition
        total_cost += np.sum((simulated_nadh_free - experimental_nadh_free) ** 2)
        total_cost += np.sum((simulated_nadh_bound - experimental_nadh_bound) ** 2)
        total_cost += np.sum((simulated_nadh_ratio - experimental_nadh_ratio) ** 2)

    return total_cost

# Initial guess for parameters
initial_params = [100, 50, 0, 0, 0, 0.01, 1, 0.01, 1, 0.01, 1]  # Adjust as needed
param_bounds = [(0, None), (0, None), ...]  # Complete with appropriate bounds

# Run the optimization
result = minimize(
    cost_function,
    initial_params,
    args=(model, tspan, data_nadhf, data_nadhb, data_bound_ratio),
    bounds=param_bounds,
    method='L-BFGS-B'
)
)

# After optimization
if result.success:
    fitted_params = result.x
    print(f"Optimized Parameters: {fitted_params}")

    filepath = '/Users/raegangiberson/Desktop/BMEG-470V_570V/Team D/optimized_parameters.json'

    try:
        with open(filepath, 'w') as file:
            json.dump(fitted_params.tolist(), file)
        print(f"Parameters saved to {filepath}")
    except Exception as e:
        print(f"Error writing to {filepath}: {e}")
else:
    print("Optimization was not successful.")
