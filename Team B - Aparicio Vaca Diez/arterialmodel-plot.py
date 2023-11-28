"""
Plots results for the paper

"""
from pylab import *
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

    # each of these is a dictionary keyed by nodes
    run1, run2, run3, run4 = data

    # applies smoothing to all values
    for run in (run1, run2, run3, run4):
        for key, values in list(run.items()):
            run[key] = smooth(values, w=10)

    #
    # Plotting Collagen
    #
    subplot(121)
    col3_1, col3_2, col3_3, col3_4 = run1['col3'], run2['col3'], run3['col3'], run4['col3']

    ps = [plot(col3_1, 'bo-'), plot(col3_2, 'ro-', label='MCL1-over-Apop'), plot(col3_3, 'b^-', label='sFas-over-Apop'),
          plot(col3_4, 'r^-', label='Collagen_3')]
    # legend( ps, ['Normal-Apop', 'MCL1-over-Apop','sFas-over-Apop','LGL-like-Apop' ], loc='best' )
    plot(col3_1, 'bo-', label='Col3_run1')
    plot(col3_2, 'ro-', label='Col3_run2')
    plot(col3_3, 'b^-', label='Col3_run3')
    plot(col3_4, 'r^-', label='Col3_run4')
    legend(loc=0)
    title(' Changes in Collagen3')
    xlabel('Time Steps')
    ylabel('Percent (%)')
    ylim((-0.1, 1.1))
    #
    # Plotting FasL and Ras
    #

    subplot(122)
    fasL1, fasL2 = run1['FasL'], run4['FasL']
    ras1, ras2 = run1['Ras'], run4['Ras']

    # ps = [ plot( fasL1, 'bo-' ), plot( fasL2, 'ro-' ), plot( ras1, 'b^-' ), plot( ras2, 'r^-' ) ]
    # legend( ps, 'Normal-FasL LGL-like-FasL Normal-Ras LGL-like-Ras'.split() , loc='lower left' )
    plot(fasL1, 'bo-', label='Normal-FasL')
    plot(fasL2, 'ro-', label='LGL-like-FasL')
    plot(ras1, 'b^-', label='Normal-Ras')
    plot(ras2, 'r^-', label='LGL-like-Ras')
    legend(loc=0)
    title(' Changes in FasL and Ras')
    xlabel('Time Steps')


if __name__ == '__main__':
    # resize this to change figure size
    figure(num=None, figsize=(14, 7), dpi=80, facecolor='w', edgecolor='k')
    make_plot()
    savefig('Figure2.png')
    show()


