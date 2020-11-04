from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

# Open the 1D grid files
f = open('HIG-19-013-Llh-1d-Grid.txt', 'r')


# Plot the grids
fig = plt.figure()
text = [[float(num) for num in line.split()] for line in f]
f.close()
text = np.array(text)
textt = np.transpose(text)
x = textt[0]
y = textt[1]
plt.plot(x,y,'.',markersize=3, color = 'blue',label="Observed")

# Find parameters
def func(X, cen2, sig1, sig2):
    z = X[:]
    return (z-cen2)**2/(sig1+sig2*(cen2-z))

def fit5para(xr, yr):

    guess = (1.4, 0.1, 0.02)

    AR, pcov = optimize.curve_fit(func, xr, yr, guess)
    return AR

ff = fit5para(x, y)
cen2,sig1,sig2 = ff[0], ff[1], ff[2]
print("\n Parameters:")
print("cen2, sig1, sig2 =", ff[0], ",", ff[1], ",", ff[2])

# Find central value and uncertainties from those parameters
print("\n Central value and uncertainties:")
def g1(y):
    return y**2 - sig2*y - sig1

a2 = fsolve(g1,1.5)
a1 = a2 + sig2
print("cen =", cen2)
print("sig2p =", a1)
print("sig2m =", a2)


x = np.arange(0,5,0.005)

y2 = (x-cen2)**2/(sig1+sig2*(-x+cen2))


plt.plot(x,y2,'-',markersize=2, color = 'g',label="Barlow's V. Gaussian Appx.")

# Choose VBF - ggH
plt.xlabel(r'$\mu_{\gamma\gamma}^{ttH}$', fontsize=20)
# plt.xlabel(r'$\mu_{ZZ}^{ggH}$', fontsize=20)
# plt.xlabel(r'$\mu_{ZZ}^{ttH}$', fontsize=20)
# plt.xlabel(r'$\mu_{ZZ}^{VH}$', fontsize=20)

plt.ylabel(r'-2 Loglikelihood', fontsize=20)
plt.title("$\mu$ from CMS-HIG-19-013 (VGaussian)")
plt.legend(loc='upper right', fontsize=12)

fig.set_tight_layout(True)


fig.savefig('mu_VBF_1D_VGaussian-fitted.pdf')
plt.show()
