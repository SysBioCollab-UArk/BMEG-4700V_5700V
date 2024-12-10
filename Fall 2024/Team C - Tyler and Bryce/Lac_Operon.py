import matplotlib.pyplot as plt
from boolean2 import Model

#
# This initial condition leads to a cycle of period 4.
# If A is set to False, a steady state is obtained.
#
#
text = """
Glucose = False
PTS = False
AC = False
cAMP = False
CAP = False
Laci = False
LacZYA = True
LacZYA_MRNA = False
LacZ = False
LacY = False
LacA = Flase
Allolactose = False
Lactose_int = False
Lactose_ex = True


PTS* = Glucose
AC* = Not PTS
cAMP* = AC
CAP* = cAMP
Laci* = Not Allolactose
LacZYA_MRNA* = CAP and lacZYA and not Laci
LacZ* = LacZYA_MRNA
LacY* = LacZYA_MRNA
LacA* = LacZYA_MRNA
Allolactose* = LacZ and Lactose_int
Lactose_int* = LacY and Lactose_ex and not PTS
"""

model = Model( text=text, mode='sync')
model.initialize()
model.iterate( steps=15 )

# the model data attribute holds the states keyed by nodes
for node in model.data:
    print(node, model.data[node])

# this is a helper function that reports the cycle lengths
# and the  index at wich the cycle started
model.report_cycles()

#
# the same thing as above but
# will not print only return the two parameters
#
print(model.detect_cycles())

#
# this is how one plots the values, delete this below
# if matplotlib is not installed
#
p1 = plt.plot( model.data["B"] , 'ob-', label='B')
p2 = plt.plot( model.data["C"] , 'sr-', label='C')
plt.legend( loc='best')
plt.ylim((-0.1,1.1))
plt.show()


