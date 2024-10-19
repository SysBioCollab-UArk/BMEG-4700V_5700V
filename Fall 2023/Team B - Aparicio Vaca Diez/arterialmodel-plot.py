"""
Plots results for the paper

"""
from matplotlib import pylab
from boolean2 import util


def smooth(data, w=0):
    "Smooths data by a moving window of width 'w'"
    fw = float(w)

    def average(index):
        return sum(data[index: index + w]) / fw

    indices = range(len(data) - w)
    out = list(map(average, indices))
    return out


def make_plot():
    # contains averaged node information based on 1000 runs
    data = util.bload('arterial_model.bin')
    print("hello data")
    # each of these is a dictionary keyed by nodes

    # applies smoothing to all values
    run = data
    print("hello run")
    for key, values in list(run.items()):
        run[key] = smooth(values, w=10)

    #
    # Plotting Collagen
    #

    col3 = run['col3']

    #ps = [plot(col3, 'col3'), label='Collagen_3')]
    # legend( ps, ['Normal-Apop', 'MCL1-over-Apop','sFas-over-Apop','LGL-like-Apop' ], loc='best' )
    plot(col3, 'col3', label='Col3')

    legend(loc=0)
    title(' Changes in Collagen3')
    xlabel('Time Steps')
    ylabel('Percent (%)')
    ylim((-0.1, 1.1))
    #
    # Plotting FasL and Ras
    #


if __name__ == '__main__':
    # resize this to change figure size
    pylab.figure(num=None, figsize=(14, 7), dpi=80, facecolor='w', edgecolor='k')
    make_plot()
    savefig('Figure2.png')
    show()


