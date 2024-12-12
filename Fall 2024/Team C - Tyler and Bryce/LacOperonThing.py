import matplotlib.pyplot as plt
from boolean2 import Model, util

text = """
Glucose = False
LactoseI = False
LactoseO = True
Allolactose = True
LacI = True
LacZ = False
LacY = False
LacZYAmRNA = True

LacZ* = LacZYAmRNA
LacY* = LacZYAmRNA
LacZYAmRNA* = not Glucose and LactoseO
Allolactose* = LacZ and LactoseI
LactoseI* = LacY and LactoseO
LacI* = not Allolactose
"""

SPECIES_TO_PLOT = ['LacZ', 'LacY', 'LacI', 'Allolactose', 'LacZYAmRNA', 'Allolactose']

##### SYNCHRONOUS UPDATING #####

model = Model(text=text, mode='sync')
model.initialize()
model.iterate(steps=15)

# the model data attribute holds the states keyed by nodes
for node in model.data:
    print(node, model.data[node])

# this is a helper function that reports the cycle lengths 
# and the  index at wich the cycle started
model.report_cycles()

# the same thing as above but
# will not print only return the two parameters
print(model.detect_cycles())    

# this is how one plots the values, delete this below if matplotlib is not installed
for sp in SPECIES_TO_PLOT:
    plt.plot(model.data[sp], marker='o', label=sp)
plt.legend( loc = "best")
plt.ylim((-0.1,1.1))
plt.title('Lac Operon: Synchronous Updating')
plt.tight_layout()

##### ASYNCHRONOUS UPDATING #####

model = Model(text, mode='async')
coll  = util.Collector()

nsims = 1000  # set the number of simulations to run to average over

for sim in range(nsims):
    model.initialize()
    model.iterate(steps=15)
    nodes = model.nodes
    # this collects states for each run
    coll.collect(states=model.states, nodes=nodes)

# averages the values for each node and  returns a dictionary keyed by nodes and a list of values, with the average
# state for each timestep
avgs = coll.get_averages(normalize=True)

plt.figure()
for sp in SPECIES_TO_PLOT:
    plt.plot(avgs[sp], marker='o', label=sp)
plt.legend( loc = "best")
plt.ylim((-0.1,1.1))
plt.title('Lac Operon: Asynchronous Updating (%d sims)' % nsims)
plt.tight_layout()

plt.show()


