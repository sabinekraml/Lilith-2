import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.optimize as optimize

# ---------- 2D Continuous Poisson Dist. with Negative Correlation ----------
# ---------- (ggF marginal - VBF dependent) ----------

## ggF
cen1 = 0.9994350282485871
## VBF
cen2 = 1.1998058252427213

# Calculate Loglikelihood

# f = open('HIGG-2017-07_Llh-2d-Grid68.txt', 'r')
f = open('HIGG-2017-07_Llh-2d-Grid95.txt', 'r')

text = [[float(num) for num in line.split()] for line in f]
f.close()
text = np.array(text)
textt = np.transpose(text)
x = textt[0]
y = textt[1]


def func(X, alpha1, alpha2, nu1, nu2, alpha):
    z1, z2 = X[:,0], X[:,1]
    A = np.exp(alpha)-1
    L2t1 = - alpha1 * (z1 - cen1) + nu1 * np.log(1 + alpha1 * (z1 - cen1) / nu1)
    L2t2a = - (alpha2 * (z2 - cen2) + nu2) * np.exp(alpha * nu1 - A * (alpha1 * (z1 - cen1) + nu1))
    L2t2b = - nu2 * np.exp((alpha - A) * nu1)
    L2t2c = nu2 * np.log(L2t2a/L2t2b)
    L2t2 = L2t2a - L2t2b + L2t2c
    L2t = -2 * (L2t1 + L2t2)
    return L2t

def fit5para(xr, yr):
    A = []
    B = []
    for i in range(0, len(xr)):
        A.append([xr[i], yr[i]])
        #B.append(2.30)
        B.append(5.99)
    A = np.array(A)
    B = np.array(B)
    guess = (2, 10, 13, 50, -0.005)
    AR, pcov = optimize.curve_fit(func, A, B, guess)
    return AR

ff = fit5para(x, y)
alpha1, alpha2, nu1, nu2, alpha = ff[0], ff[1], ff[2], ff[3], ff[4]

x = np.exp(alpha)-1
corr = nu1*nu2*x/np.sqrt(nu1*nu2*(1+nu2*(np.exp(nu1*x**2)-1)))

gam1 = alpha1/nu1
gam2 = alpha2/nu2

print("\n Central value and uncertainties:")

def g1(sig1p):
    return gam1*sig1p - np.log(1 + gam1*sig1p) - 1/(2*nu1)

def g2(sig1m):
    return np.exp(-gam1*(sig1m + sig1p)) - (1 - gam1*sig1m)/(1 + gam1*sig1p)

sig1p = fsolve(g1,1)
sig1m = fsolve(g2,1)

print("cen1, sig1p, sig1m =", cen1,sig1p[0], sig1m[0])

def s1(sig2p):
    return gam2*sig2p - np.log(1 + gam2*sig2p) - 1/(2*nu2)

def s2(sig2m):
    return np.exp(-gam2*(sig2m + sig2p)) - (1 - gam2*sig2m)/(1 + gam2*sig2p)

sig2p = fsolve(s1,1)
sig2m = fsolve(s2,1)

print("cen2, sig2p, sig2m =", cen2,sig2p[0], sig2m[0])

print("correlation =", corr)

A = np.exp(alpha)-1

def Loglikelihood(z1,z2):
    L2t1 = - alpha1 * (z1 - cen1) + nu1 * np.log(1 + alpha1 * (z1 - cen1) / nu1)
    L2t2a = - alpha2 * (z2 - cen2 + 1 / gam2) * np.exp(alpha * nu1 - A * alpha1 * (z1 - cen1 + 1/gam1))
    L2t2b = - alpha2 * ( 1 / gam2) * np.exp(alpha * nu1 - A * alpha1 * ( 1 / gam1))
    L2t2c = nu2 * np.log(L2t2a/L2t2b)
    L2t2 = L2t2a - L2t2b + L2t2c
    L2t = -2 * (L2t1 + L2t2)
    return L2t

# Grid calculations
xmax, ymax = 3.3, 3.3
(a, b) = (0.01, 0.01)
x2 = np.arange(-0.5, xmax, a, dtype=np.float)
y2 = np.arange(-0.5, ymax, b, dtype=np.float)
X2, Y2 = np.meshgrid(x2, y2)
Z = Loglikelihood(X2,Y2)
print("Minimum LogLikelihood =",Z.min())

# Open the 2D grid files
f = open('HIGG-2017-07_Llh-2d-GridFull.txt', 'r')

# Plot the grids
text = [[float(num) for num in line.split()] for line in f]
f.close()
text = np.array(text)
textt = np.transpose(text)
expx = textt[0]
expy = textt[1]


######################################################################
# Plot routine
######################################################################

print "***** plotting *****"

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 15
matplotlib.rcParams['ytick.major.pad'] = 15

fig = plt.figure()
ax = fig.add_subplot(111)

plt.minorticks_on()
plt.tick_params(labelsize=20, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

# Plotting the 68%, 95% CL regions

ax.contourf(X2,Y2,Z,[1e-10,2.3,5.99],colors=['gray','silver','#ffff00'], \
              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

# best fit point
plt.plot([1.0],[1.2], '*', c='w', ms=10)

# official ATLAS contour
plt.plot(expx,expy, '.', c='b', label='ATLAS official')
plt.legend(loc='lower right', fontsize=12)

# SM value
#plt.plot([1],[1], '+', c='k', ms=10)

# Title, labels, color bar...
ax.set_aspect(0.7)
plt.title("               Lilith-2.0, DB 19.06", fontsize=14.5, ha="left")
plt.xlim(-1,4)
plt.ylim(-1,4)
plt.xlabel(r'$\mu({\rm ggH},\tau\tau)$',fontsize=25)
plt.ylabel(r'$\mu({\rm VBF},\tau\tau)$',fontsize=25)
plt.text(-0.7, 3.45, 'Data from ATLAS-HIGG-2017-07:', fontsize=13)
plt.text(-0.7, 3.17, '2D Poisson fitted to ATLAS 95% CL contour', fontsize=13)

fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig("HIGG-2017-07-mu-2d-Poisson.pdf")


