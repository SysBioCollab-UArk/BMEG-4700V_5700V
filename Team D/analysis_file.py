import pandas as pd
import numpy as np
from pysb.simulator import ScipyOdeSimulator
from scipy.optimize import minimize
from redox_ratio import model  # Make sure this is the correct path to your model

# Load experimental data
exp_data_fig2 = pd.read_excel('/mnt/data/elife-73808-fig2-data1-v2.xlsx')
# ... load other data files similarly

# Define time span for the simulation
tspan = np.linspace(0, len(exp_data_fig2) - 1, len(exp_data_fig2))


# Define a function to simulate the model
def simulate_model(params, model, tspan):
    # Update model parameters
    for i, param in enumerate(model.parameters):
        if i < len(params):
            param.value = params[i]

    # Simulate the model
    sim = ScipyOdeSimulator(model, tspan)
    return sim.run()


# Define the cost function
def cost_function(params, model, experimental_data, tspan):
    simulation_results = simulate_model(params, model, tspan)
    simulated_data = simulation_results.observables['NADH_free']

    cost = 0
    for data_df in experimental_data:
        # Replace 'NADHf' and 'NADHb' with the actual column names in your data
        experimental_data_f = data_df['NADHf'].values
        experimental_data_b = data_df['NADHb'].values
        cost += np.sum((simulated_data - experimental_data_f) ** 2)
        cost += np.sum((simulated_data - experimental_data_b) ** 2)

    return cost


# Combine experimental data into a list
experimental_data_list = [exp_data_fig2]  # Add other data frames to this list

# Initial guess for parameters
initial_params = [100, 50, 0, 0, 0, 0.01, 1, 0.01, 1, 0.01, 1]  # Adjust as needed
param_bounds = [(0, None), (0, None), ...]  # Define bounds

# Run the optimization
result = minimize(
    cost_function,
    initial_params,
    args=(model, experimental_data_list, tspan),
    bounds=param_bounds,
    method='L-BFGS-B'
)

# Evaluate the optimized model
if result.success:
    fitted_params = result.x
    optimized_simulation = simulate_model(fitted_params, model, tspan)

    # Calculate the discrepancy
    discrepancy = cost_function(fitted_params, model, experimental_data_list, tspan)
    print("Discrepancy between model and experimental data:", discrepancy)
else:
    raise ValueError(result.message)
