from pysb import *
from pysb.simulator import ScipyOdeSimulator
from scipy.optimize import minimize
import numpy as np
import pandas as pd

# Assuming 'model' is your PySB model defined in redox_ratio.py
from redox_ratio import model

# Load your experimental data into pandas DataFrames
experimental_data_fig2 = pd.read_excel('elife-73808-fig2-data1-v2.xlsx')
experimental_data_fig5 = pd.read_excel('elife-73808-fig5-data1-v2.xlsx')
experimental_data_fig6 = pd.read_excel('elife-73808-fig6-data1-v2.xlsx')
experimental_data_fig7 = pd.read_excel('elife-73808-fig7-data1-v2.xlsx')
experimental_data_fig8 = pd.read_excel('elife-73808-fig8-data1-v2.xlsx')

# Combine the data into a list for easy iteration in the cost function
experimental_data_list = [
    experimental_data_fig2,
    experimental_data_fig5,
    experimental_data_fig6,
    experimental_data_fig7,
    experimental_data_fig8
]

# Define the time span for the simulation based on your experimental data
tspan = np.linspace(0, len(experimental_data_fig2) - 1,
                    len(experimental_data_fig2))  # Adjust according to the time points in your data


# Define the cost function that will be minimized
def cost_function(params, model, experimental_data_list, tspan):
    # Update model parameters
    for i, param in enumerate(model.parameters):
        param.value = params[i]

    # Simulate the model
    sim = ScipyOdeSimulator(model, tspan)
    simulation_results = sim.run()

    # Initialize cost
    cost = 0

    # Compare each observable with its corresponding experimental data
    for data_df in experimental_data_list:
        for observable in model.observables:
            # Assuming the observable name matches the column name in your data_df
            if observable.name in data_df.columns:
                simulated_data = simulation_results.observables[observable.name]
                experimental_data = data_df[observable.name].values
                # Update the cost based on the difference between simulation and experiment
                cost += np.sum((simulated_data - experimental_data) ** 2)

    return cost


# Initial guess for parameters
# Replace with your model's initial parameter values
initial_params = [100, 50, 0, 0, 0, 0.01, 1, 0.01, 1, 0.01, 1]

# Bounds for parameters
# Replace with appropriate bounds or set to None if unknown
param_bounds = [(0, None), (0, None), ...]  # Add the rest of the bounds for all parameters

# Run the optimization
result = minimize(
    cost_function,
    initial_params,
    args=(model, experimental_data_list, tspan),
    bounds=param_bounds,
    method='L-BFGS-B'  # Example method, choose the one that suits your problem
)

# Check the result
if result.success:
    fitted_params = result.x
    print("Fitted parameters:", fitted_params)
else:
    raise ValueError(result.message)
