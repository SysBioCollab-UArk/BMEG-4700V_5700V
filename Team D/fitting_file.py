import pandas as pd
import numpy as np
from pysb.simulator import ScipyOdeSimulator
from scipy.optimize import minimize
from redox_ratio import model

# Load experimental data into pandas DataFrames
experimental_data_fig2 = pd.read_excel('/mnt/data/elife-73808-fig2-data1-v2.xlsx')  # Example: oxamate condition
experimental_data_fig5 = pd.read_excel('/mnt/data/elife-73808-fig5-data1-v2.xlsx')  # Example: rot condition
experimental_data_fig6 = pd.read_excel('/mnt/data/elife-73808-fig6-data1-v2.xlsx')  # Additional condition
experimental_data_fig7 = pd.read_excel('/mnt/data/elife-73808-fig7-data1-v2.xlsx')  # Additional condition
experimental_data_fig8 = pd.read_excel('/mnt/data/elife-73808-fig8-data1-v2.xlsx')  # Additional condition

# Combine the data into a list for easy iteration in the cost function
experimental_data_list = [
    experimental_data_fig2,
    experimental_data_fig5,
    experimental_data_fig6,
    experimental_data_fig7,
    experimental_data_fig8
]

# Define the time span for the simulation based on your experimental data
tspan = np.linspace(0, len(experimental_data_fig2) - 1, len(experimental_data_fig2))

# Define the cost function that will be minimized
def cost_function(params, model, experimental_data_list, tspan):
    total_cost = 0

    # Update model parameters
    for i, param in enumerate(model.parameters):
        if i < len(params):
            param.value = params[i]

    # Iterate over each data set
    for data_df in experimental_data_list:
        # Update model for current condition if necessary
        # ...

        # Simulate the model
        sim = ScipyOdeSimulator(model, tspan)
        simulation_results = sim.run()

        # Compare simulation with experimental data
        simulated_nadh_free = simulation_results.observables['NADH_free']
        simulated_nadh_bound = simulation_results.observables['NADH_bound']
        simulated_nadh_ratio = simulation_results.observables['NADH_ratio']

        experimental_nadh_free = data_df['NADHf'].values  # Replace 'NADHf' with actual column name
        experimental_nadh_bound = data_df['NADHb'].values  # Replace 'NADHb' with actual column name
        experimental_nadh_ratio = data_df['NADH_ratio'].values  # Replace 'NADH_ratio' with actual column name

        total_cost += np.sum((simulated_nadh_free - experimental_nadh_free) ** 2)
        total_cost += np.sum((simulated_nadh_bound - experimental_nadh_bound) ** 2)
        total_cost += np.sum((simulated_nadh_ratio - experimental_nadh_ratio) ** 2)

    return total_cost

# Initial guess for parameters
initial_params = [100, 50, 0, 0, 0, 0.01, 1, 0.01, 1, 0.5]  # Adjust as needed
param_bounds = [(0, None), (0, None), ...]  # Complete with appropriate bounds

# Run the optimization
result = minimize(
    cost_function,
    initial_params,
    args=(model, experimental_data_list, tspan),
    bounds=param_bounds,
    method='L-BFGS-B'
)

# Check the result
if result.success:
    fitted_params = result.x
    print("Fitted parameters:", fitted_params)
else:
    raise ValueError(result.message)
