from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

# Open the 1D grid files
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
cen1, sig1m, sig1p = 1.04, 0.20, 0.245
#
# Poisson method
# Calculate parameters
def f(t):
    return np.exp(-t*(sig1m + sig1p)) - (1 - t*sig1m)/(1 + t*sig1p)

gam = fsolve(f,1.001)
L = 1/(2*(gam*sig1p - np.log(1 + gam*sig1p)))
alpha = L*gam
print(alpha,L)
# Range of validation
n = 5
x2 = np.arange(1.8,n,0.02)
# Log Likelihood
y2 = -2*(-alpha*(x2 - cen1) + L*np.log(1 + alpha*(x2 - cen1)/L))
# # Plot
plt.plot(x2,y2,'.',markersize=2, color = 'g',label="Barlow's Poisson Appx.")

# # VGaussian method
# # Calculate parameters
# sig1 = (2*sig1m*sig1p)/(sig1m+sig1p)
# sig2 = (sig1p-sig1m)/(sig1p+sig1m)
# # Range of validation
# # x2 = np.arange(0,n,0.01)
# # Log Likelihood
# y2 = (x2-cen1)**2/(sig1+sig2*(x2-cen1))**2

f = open('CONF-2020-045-Llh-1d-Grid-added.txt', 'w')
for i in range(len(x2)):
    f.write(str(np.round((x2[i]),7))+'\t'+str(y2[i])+'\n')
f.close()

# Plot
# plt.plot(x2,y2,'-',markersize=2, color = 'r',label="Barlow's V. Gaussian Appx.")

# # Gaussian method
# # Range of validation
# x2 = np.arange(0,n,0.005)
# # Log Likelihood
# # Left side
# x3 = np.arange(0,cen1,0.005)
# y3 = (x2 - cen1) ** 2 / sig1m ** 2
# # Right side
# x4 = np.arange(cen1,2,0.005)
# y4 = (x2 - cen1) ** 2 / sig1p ** 2
# y2 = y3 + y4

# Plot
# plt.plot(x2,y2,'-',markersize=2, color = 'y',label="Gaussian Appx.")

plt.xlabel(r'$\mu_{ZZ}^{ggH}$', fontsize=20)
plt.ylabel(r'-2 Loglikelihood', fontsize=20)
plt.title("$\mu$ from ATLAS-CONF-2020-045")
plt.legend(loc='upper right', fontsize=12)

fig.set_tight_layout(True)

# Choose VBF - ggH
fig.savefig('Grids_added_Poisson_method.pdf')

plt.show()
