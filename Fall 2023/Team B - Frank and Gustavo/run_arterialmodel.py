"""
LGL simulator

It is also a demonstration on how the collector works

"""
import boolean2
import os
from boolean2 import Model, util
from random import choice
from matplotlib import pylab
import matplotlib.pyplot as plt

# ocasionally randomized nodes
TARGETS = set("Col1".split())


def new_getvalue(state, name, p):
    """
    Called every time a node value is used in an expression.
    It will override the value for the current step only.
    Returns random values for the node states
    """
    global TARGETS
    value = util.default_get_value(state, name, p)

    if name in TARGETS:
        # pick at random from True, False and original value
        return choice([True, False, value])
    else:
        return value


def run(text, nodes, repeat, steps):
    """
    Runs the simulation and collects the nodes into a collector,
    a convenience class that can average the values that it collects.
    """
    coll = util.Collector()

    for i in range(repeat):
        engine = Model(mode='async', text=text)
        #engine.RULE_GETVALUE = new_getvalue
        # minimalist initial conditions, missing nodes set to false
        engine.initialize()  # missing=util.false)
        engine.iterate(steps=steps)
        coll.collect(states=engine.states, nodes=nodes)
        print('repeat no. i=%s' % i)

    print('- completed')
    avgs = coll.get_averages(normalize=True)

    return avgs


if __name__ == '__main__':
    # read in the text
    #text = file('temp.txt').read()


    # Get the absolute path of the directory containing the script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the absolute path to the 'temp.txt' file
    file_path = os.path.join(script_dir, 'temp.txt')

    # Read the contents of the file
    try:
        with open(file_path) as f:
            text = f.read()
            print("File Content: {}".format(repr(text)))
    except IOError as e:
        print("Error: {}".format(e))



    #print("Current Working Directory:", os.getcwd())


    # the nodes of interest that are collected over the run
    # NODES  = 'Apoptosis STAT3 FasL Ras'.split()

    # this collects the state of all nodes
    NODES = boolean2.all_nodes(text)
    print()
    #
    # raise this for better curves (will take about 2 seconds per repeat)
    # plots were made for REPEAT = 1000, STEPS=150
    #
    REPEAT = 20
    STEPS = 50
    data = []
    print('- starting simulation with REPEAT=%s, STEPS=%s' % (REPEAT, STEPS))
    avgs = run(text=text, repeat=REPEAT, nodes=NODES, steps=STEPS)
    data.append(avgs)
    fname = 'arterial_model.bin'
    util.bsave(data, fname=fname)



    col3 = avgs['Col3']
    col1 = avgs['Col1']
    print('extracted variables for graph')
    print(col3[1])
    print(col3[50])


    plt.plot(col3,label='Col 3')
    print('plotting col3')
    plt.plot(col1, label='Col 1')
    print('plotting col1')

    # Add legend
    plt.legend()

    # Add labels to the x and y axes
    plt.xlabel('Round')
    plt.ylabel('Count')

    # Add a title to the plot
    plt.title('Col 1 and 3')

    # Show the plot
    plt.show()

    print('finished graph')

