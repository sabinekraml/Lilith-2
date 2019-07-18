from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import sys

# Open the 1D grid files

# Choose VBF - ggH - ttH - VH
# f = open('HIGG-2017-07_Llh-2d-Grid68.txt', 'r')
f = open('HIGG-2017-07_Llh-2d-Grid95.txt', 'r')

# Plot the grids
fig = plt.figure()
text = [[float(num) for num in line.split()] for line in f]
f.close()
text = np.array(text)
textt = np.transpose(text)
x = textt[0]
y = textt[1]

# Find parameters
def func(X, sig1p, sig1m, sig2p, sig2m, p):
    z1, z2 = X[:,0], X[:,1]
    z10, z20 = 0.9994350282485871,	1.1998058252427213
    V1 = sig1p * sig1m
    V1e = sig1p - sig1m
    V2 = sig2p * sig2m
    V2e = sig2p - sig2m
    V1f = V1 + V1e * (z1 - z10)
    V2f = V2 + V2e * (z2 - z20)
    L2t = 1 / (1 - p ** 2) * (
            (z1 - z10) ** 2 / V1f - 2 * p * (z1 - z10) * (z2 - z20) / np.sqrt(V1f * V2f) + (z2 - z20) ** 2 / V2f)
    return L2t

def fit5para(xr, yr):
    A = []
    B = []
    for i in range(0, len(xr)):
        A.append([xr[i], yr[i]])
        # B.append(2.30)
        B.append(5.99)
    A = np.array(A)
    B = np.array(B)
    guess = (1, 1, 1, 1, -0.5)

    AR, pcov = optimize.curve_fit(func, A, B, guess)
    return AR

ff = fit5para(x, y)
sig1p, sig1m, sig2p, sig2m, p = ff[0], ff[1], ff[2], ff[3], ff[4]
print("\n Parameters:")
print("sig1p, sig1m, sig2p, sig2m, corr =", ff[0], ",", ff[1], ",", ff[2], ",", ff[3], ",", ff[4])

#sys.exit("end here")

# Open the 1D grid files
f = open('HIGG-2017-07_Llh-2d-GridFull.txt', 'r')

# Plot the grids
text = [[float(num) for num in line.split()] for line in f]
f.close()
text = np.array(text)
textt = np.transpose(text)
x = textt[0]
y = textt[1]
plt.plot(x,y,'.',markersize=3, color = 'red', label="Observed")

# Loglikelihood Calculation
def Loglikelihood(z1,z2):
    z10, z20 = 0.9994350282485871, 1.1998058252427213
    V1 = sig1p*sig1m
    V1e = sig1p - sig1m
    V2 = sig2p * sig2m
    V2e = sig2p - sig2m
    V1f = V1 + V1e*(z1-z10)
    V2f = V2 + V2e * (z2 - z20)
    L2t = 1/(1-p**2)*((z1-z10)**2/V1f-2*p*(z1-z10)*(z2-z20)/np.sqrt(V1f*V2f)+(z2-z20)**2/V2f)
    return L2t

# Grid calculations
xmax, ymax = 3.3, 3.3
(a, b) = (0.01, 0.01)
x2 = np.arange(-0.5, xmax, a, dtype=np.float)
y2 = np.arange(-0.5, ymax, b, dtype=np.float)
X2, Y2 = np.meshgrid(x2, y2)
Z = Loglikelihood(X2,Y2)
print("Minimum Log likelihood =",Z.min())
plt.contour(X2,Y2,Z,[0.00001,2.30,5.99])


# Print data
g = open('2D_CL_VGaussian.txt', 'w')
for i in range (len(x2)):
    for j in range (len(y2)):
        g.write(str(x2[i]) + " " + str(y2[j]) + " " + str(Z[j][i]) + '\n')
g.close


# Plot
plt.legend(loc='upper right', fontsize=12)
plt.xlabel(r'$\mu_{ZZ}^{VBF}$', fontsize=30)
plt.ylabel(r'$\mu_{ZZ}^{ggF}$', fontsize=30)
# plt.title("$\mu$ from ATLAS-HIGG-2017-07 (VGaussian, 68% CL, fitted corr)")
plt.title("$\mu$ from ATLAS-HIGG-2017-07 (VGaussian, 95% CL, fitted corr)")
fig.set_tight_layout(True)

# fig.savefig('mu_ggH-VBF_2D_VGaussian-fitted-corr-68CL.pdf')
# fig.savefig('mu_ggH-VBF_2D_VGaussian-fitted-corr-95CL.pdf')
