import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import STU_2HDM as calc
from iminuit import Minuit
from iminuit.cost import LeastSquares

# our line model, unicode parameter names are supported :)
def line(x, a, b):
    return a+x*b

# generate random toy data with random offsets in y
np.random.seed(1)
data_x = np.linspace(0, 1, 10)
data_yerr = 0.1  # could also be an array
data_y = line(data_x, 1, 2) + data_yerr * np.random.randn(len(data_x))


# iminuit contains a LeastSquares class to conveniently generate a least-squares cost function.
# We will revisit how to write this by hand in a later section.
least_squares = LeastSquares(data_x, data_y, data_yerr, line)

m = Minuit(least_squares, a=0, b=0)  # starting values for α and β

print("start migrad")
m.migrad()  # finds minimum of least_squares function
#m.hesse()   # accurately computes uncertainties
print("start minos")
m.minos()

## draw data and fitted line
#plt.errorbar(data_x, data_y, data_yerr, fmt="o", label="data")
#plt.plot(data_x, line(data_x, *m.values), label="fit")

# draw three contours with 68%, 90%, 99% confidence level
#m.draw_mncontour("a", "b", cl=(0.68, 0.9, 0.99));
print("start mncontour")
test = m.mncontour("a","b")
#print("test = ", test)
#plt.show()
