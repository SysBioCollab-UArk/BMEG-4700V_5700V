# BMEG-470V_570V

import matplotlib.pyplot as plt
import boolean2
from boolean2 import util

# def convert_dict_to_text(dictionary):
#     text = ""
#     for key, value in dictionary.items():
#         text += f"{key} = {value}\n"
#     return text.strip()


def load_rules_from_file(file_path):
    initial_state = {}
    rules = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                if line[0] == '1':  # rules
                    variable, rule = map(str.strip, line.split('='))
                    rules[variable] = rule
                else:  # initial states
                    variable, state = map(str.strip, line.split('='))
                    initial_state[variable] = state
    return initial_state, rules


# initial_state = {
#         'EGF': True,
#         'DUSP': False,
#         'Cyt_Ca': False,
#         'RAS': False,
#         'RAF': False,
#         'MEK': False,
#         'ERK': False,
#         'TIMP': False,
#         'MMP': False,
#         'Collagen': False,
#         'PI3K': False,
#         'AKT': False,
#         'mTOR': False,
#         'BNP': False
#     }

# initial_state_text = convert_dict_to_text(initial_state)


def simulate_boolean_model(initial_state, rules, num_steps, num_sims):
    # variable_values = {variable: [] for variable in initial_state}
    #
    # return variable_values

    # num_steps = 10

    coll = util.Collector()

    model_text = ''
    for item in initial_state.items():
        model_text += '%s = %s\n' % (item[0], str(item[1]))
    model_text += '\n'
    for item in rules.items():
        model_text += '%s = %s\n' % (item[0], item[1])

    for n in range(num_sims):
        model = boolean2.Model(text=model_text, mode='async')
        model.initialize()
        print('sim %d' % n)
        model.iterate(steps= num_steps)
        # print(type(model.states))
        # print(len(model.states))
        # for i in range(len(model.states)):
        #     print(model.states[i])

        coll.collect(states= model.states, nodes= model.nodes)

    avgs = coll.get_averages(normalize=True)
    return avgs

def plot_simulation_results(variable_values, num_steps):

    plt.figure(figsize=(16,8))
    variables_dict = {key.lower(): key for key in variable_values.keys()}

    for variable in sorted([key.lower() for key in variable_values.keys()]):
        fmt= '-'
        ms = 8

        if variables_dict[variable] == 'Collagen' or variables_dict[variable] == 'BNP':
            fmt='o-'
        elif variables_dict[variable]== 'Cyt_Ca':
            fmt='s-'
        elif variables_dict[variable]== 'EGF':
            fmt= '^-'
        elif variables_dict[variable]== 'DUSP':
            fmt = 'v-'
            ms = 12
        plt.plot([i for i in range(num_steps + 1)], variable_values[variables_dict[variable]], fmt, ms = ms,lw=2,
                 label=variables_dict[variable])


    plt.xlabel('Step', fontsize= 20)
    plt.ylabel('Value',fontsize= 20)
    plt.xticks(fontsize= 16)
    plt.yticks(fontsize= 16)
    plt.legend(loc='best',fontsize= 16)
    plt.title('Boolean Model Simulation',fontsize= 20)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    rules_file = 'Class_Project.txt'
    num_sims = 1000
    num_steps = 10
    initial_state, rules = load_rules_from_file(rules_file)

    # print(initial_state)
    # print(rules)
    # print(num_steps)
    simulation_results = simulate_boolean_model(initial_state, rules, num_steps, num_sims)
    # print(simulation_results)
    plot_simulation_results(simulation_results, num_steps)
