"""
LGL simulator

It is also a demonstration on how the collector works

"""
import boolean2
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

    print('- completed')
    avgs = coll.get_averages(normalize=True)

    return avgs


if __name__ == '__main__':
    # read in the text
    text = file('temp.txt').read()

    # the nodes of interest that are collected over the run
    # NODES  = 'Apoptosis STAT3 FasL Ras'.split()

    # this collects the state of all nodes
    NODES = boolean2.all_nodes(text)
    print()
    #
    # raise this for better curves (will take about 2 seconds per repeat)
    # plots were made for REPEAT = 1000, STEPS=150
    #
    REPEAT = 10
    STEPS = 50
    data = []
    print('- starting simulation with REPEAT=%s, STEPS=%s' % (REPEAT, STEPS))
    avgs = run(text=text, repeat=REPEAT, nodes=NODES, steps=STEPS)
    data.append(avgs)
    fname = 'arterial_model.bin'
    util.bsave(data, fname=fname)



    col3 = avgs['Col3']
    col1 = avgs['Col1']
    variable_type = type(col3)
    print("The type of avgs is:", variable_type)
    print("The value using square bracket notation", col3[3])

    plt.plot(col3,label='Col 3')
    plt.plot(col1, label='Col 1')

    # Add legend
    plt.legend()

    # Add labels to the x and y axes
    plt.xlabel('Round')
    plt.ylabel('Count')

    # Add a title to the plot
    plt.title('Col 1 and 3')

    # Show the plot
    plt.show()


