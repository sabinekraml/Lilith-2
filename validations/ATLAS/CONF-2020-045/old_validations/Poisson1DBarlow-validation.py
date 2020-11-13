from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

# Open the 1D grid files

# Choose VBF - ggH
# f = open('HIGG-2016-22_VBF-1d-Grid.txt', 'r')
f = open('CONF-2020-045-Llh-1d-Grid.txt', 'r')

fig = plt.figure()
# Plot the grids
text = [[float(num) for num in line.split()] for line in f]
f.close()
text = np.array(text)
textt = np.transpose(text)
x = textt[0]
y = textt[1]
plt.plot(x,y,'.',markersize=3, color = 'blue',label="Observed")

#Input data

# #Data from Table 9
# #ggF
# cen1, sig1m, sig1p = 1.11, np.sqrt(0.2**2 + 0.06**2 + 0.04**2), np.sqrt(0.22**2+0.07**2+0.04**2)
# #VBF
# cen2, sig2m, sig2p = 4.0, np.sqrt(1.4**2 + 2*0.3**2), np.sqrt(1.7**2 + 2*0.3**2)

#Slightly Fixed the data from Table 9
#ggF
# cen1, sig1m, sig1p = 1.11+0.01, np.sqrt(0.2**2 + 0.06**2 + 0.04**2)+0.01, np.sqrt(0.22**2+0.07**2+0.04**2)+0.01
#VBF
cen1, sig1m, sig1p = 1.04, 0.20, 0.24


def f(t):
    # Choose VBF - ggH
    # return np.exp(-t*(sig2m + sig2p)) - (1 - t*sig2m)/(1 + t*sig2p)
    return np.exp(-t*(sig1m + sig1p)) - (1 - t*sig1m)/(1 + t*sig1p)

A = fsolve(f,1.001)
print(A)

gam = A

# Choose VBF - ggH
# L = 1/(2*(gam*sig2p - np.log(1 + gam*sig2p)))
L = 1/(2*(gam*sig1p - np.log(1 + gam*sig1p)))
alpha = L*gam

print(alpha,L)
# Choose VBF - ggH
# x = np.arange(0.2,8.1,0.005)
# y2 = -2*(-alpha*(x - cen2) + L*np.log(1 + alpha*(x - cen2)/L))
x = np.arange(0,2,0.005)
y2 = -2*(-alpha*(x - cen1) + L*np.log(1 + alpha*(x - cen1)/L))


plt.plot(x,y2,'-',markersize=2, color = 'g',label="Barlow's Poisson Appx.")

# Choose VBF - ggH
# plt.xlabel(r'$\mu_{ZZ}^{VBF}$', fontsize=20)
plt.xlabel(r'$\mu_{ZZ}^{ggH}$', fontsize=20)


plt.ylabel(r'-2 Loglikelihood', fontsize=20)
plt.title("$\mu$ from ATLAS-CONF-2020-045 (Poisson)")
plt.legend(loc='upper right', fontsize=12)

fig.set_tight_layout(True)

# Choose VBF - ggH
fig.savefig('mu_VBF_1D_Poisson-validation.pdf')
# fig.savefig('mu_ggF_1D_Poisson-fixed.pdf')
plt.show()