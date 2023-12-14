import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
from pysb.simulator import ScipyOdeSimulator
from redox_ratio import model
from sklearn.metrics import mean_squared_error  # Import mean_squared_error

# Load optimized parameters from the specified file and convert to NumPy array
with open('/Users/raegangiberson/Desktop/BMEG-470V_570V/Team D/optimized_parameters.json', 'r') as file:
    fitted_params = np.array(json.load(file))  # Convert to NumPy array

# Load experimental data from different sheets
data_nadhf = pd.read_excel('/mnt/data/elife-73808-fig2-data1-v2.xlsx', sheet_name='nadhf')
data_nadhb = pd.read_excel('/mnt/data/elife-73808-fig2-data1-v2.xlsx', sheet_name='nadhb')
data_bound_ratio = pd.read_excel('/mnt/data/elife-73808-fig2-data1-v2.xlsx', sheet_name='bound ratio')

# Define time span for the simulation
tspan = np.linspace(0, len(data_nadhf) - 1, len(data_nadhf))

# Function to simulate the model
def simulate_model(params):
    for i, param in enumerate(model.parameters):
        if i < len(params):
            param.value = params[i]
    sim = ScipyOdeSimulator(model, tspan)
    return sim.run()

# Run the simulation with the optimized parameters
simulation_results = simulate_model(fitted_params)

# Analyze and plot results for each condition
conditions = ['oxamate', 'rot', 'oligo', 'fccp']
for condition in conditions:
    # Calculate RMSE for NADH free and bound
    rmse_free = np.sqrt(mean_squared_error(data_nadhf[condition], simulation_results.observables['NADH_free'].all))  # Use .all
    rmse_bound = np.sqrt(mean_squared_error(data_nadhb[condition], simulation_results.observables['NADH_bound'].all))  # Use .all

    # Print RMSE results
    print(f"RMSE for NADH Free ({condition}): {rmse_free}")
    print(f"RMSE for NADH Bound ({condition}): {rmse_bound}")

    # Plotting results
    plt.figure(figsize=(15, 5))

    # NADH Free
    plt.subplot(1, 3, 1)
    plt.plot(tspan, simulation_results.observables['NADH_free'].all, label='Model - NADH Free')  # Use .all
    plt.scatter(range(len(data_nadhf)), data_nadhf[condition], color='red', label='Experimental - NADH Free')
    plt.title(f'NADH Free - {condition}')
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.legend()

    # NADH Bound
    plt.subplot(1, 3, 2)
    plt.plot(tspan, simulation_results.observables['NADH_bound'].all, label='Model - NADH Bound')  # Use .all
    plt.scatter(range(len(data_nadhb)), data_nadhb[condition], color='red', label='Experimental - NADH Bound')
    plt.title(f'NADH Bound - {condition}')
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.legend()

    # NADH Ratio
    plt.subplot(1, 3, 3)
    nadh_ratio = simulation_results.observables['NADH_free'].all / simulation_results.observables['NADH_bound'].all  # Use .all
    plt.plot(tspan, nadh_ratio, label='Model - NADH Ratio')
    plt.scatter(range(len(data_bound_ratio)), data_bound_ratio[condition], color='red', label='Experimental - NADH Ratio')
    plt.title(f'NADH Ratio - {condition}')
    plt.xlabel('Time')
    plt.ylabel('Ratio')
    plt.legend()

    plt.tight_layout()
    plt.show()
