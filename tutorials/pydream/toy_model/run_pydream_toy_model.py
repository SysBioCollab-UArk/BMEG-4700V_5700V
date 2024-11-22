from tutorials.pydream.toy_model.toy_model import model
from param_calibration import *
from pysb.simulator import ScipyOdeSimulator

solver = ScipyOdeSimulator(model)
sim_protocol = SimulationProtocol(solver)

custom_priors = {'k_5': ('uniform', 2)}
no_sample = ['k_6']
obs_labels = {'complex_2': 'complex 2', 'product_total': 'product'}

exp_data_file = 'toy_model_synth_data.csv'

if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      sim_protocol,
                                      priors=custom_priors,
                                      no_sample=no_sample)

    calibrator.run(niterations=5000, nchains=5, obs_labels=obs_labels, plot_results=True,
                   plot_tc_args={'separate_plots': False})
