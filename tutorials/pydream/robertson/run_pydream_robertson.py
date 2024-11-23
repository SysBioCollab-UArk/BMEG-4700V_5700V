from tutorials.pydream.robertson.robertson import model
from param_calibration import *
from pysb.simulator import ScipyOdeSimulator

solver = ScipyOdeSimulator(model)
sim_protocol = SimulationProtocol(solver)

custom_priors = {'k1': ('uniform', 6), 'k2': ('uniform', 6), 'k3': ('uniform', 6)}
no_sample = ['A_0', 'B_0', 'C_0']

exp_data_file = 'robertson_synth_data.csv'

if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      sim_protocol,
                                      priors=custom_priors,
                                      no_sample=no_sample)

    calibrator.run(niterations=50000, nchains=5, plot_results=True)
