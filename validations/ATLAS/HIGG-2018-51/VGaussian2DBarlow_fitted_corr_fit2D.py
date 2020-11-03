from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize as optimize
import sys

# Open the 1D grid files

# Choose VBF - ggH - ttH - VH
# f = open('HIGG_2016-22_Llh-2d-Grid68.txt', 'r')
f = open('HIGG-2018-51-Llh-2d-Grid95.txt', 'r')

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
    z10, z20 = 0.954166667,	1.079328859
    # z10, z20 = 0.809361060171, 2.03669874783
    # # case 3
    # z10, z20 = 0.809358651978, 2.03723510311
    # z10, z20 = 0.81, 2
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
    guess = (0.6, 0.5, 0.6, 0.5, 0.3)

    AR, pcov = optimize.curve_fit(func, A, B, guess)
    return AR

ff = fit5para(x, y)
sig1p, sig1m, sig2p, sig2m, p = ff[0], ff[1], ff[2], ff[3], ff[4]

print(p)
print("\n Parameters:")
print("sig1p, sig1m, sig2p, sig2m, corr =", ff[0], ",", ff[1], ",", ff[2], ",", ff[3], ",", ff[4])
# sig1p, sig1m, sig2p, sig2m, p = 0.186957153857, 0.176141692478, 0.606157414018, 0.526259988038, -0.27
# sig1p, sig1m, sig2p, sig2m, p = 0.19777365, 0.18695455, 0.68383827, 0.60539335, -0.27

# # ggF
# cen1, sig1m, sig1p = 0.81, 0.18, 0.19
# # # VBF
# cen2, sig2m, sig2p = 2, 0.5, 0.6
# p = -0.27

#sys.exit("end here")

# # Open the 1D grid files
# f = open('HIGG_2016-21_Llh-2d-Grid.txt', 'r')
#
# # Plot the grids
# text = [[float(num) for num in line.split()] for line in f]
# f.close()
# text = np.array(text)
# textt = np.transpose(text)
# x = textt[0]
# y = textt[1]
# plt.plot(x,y,'.',markersize=3, color = 'red', label="Observed")

# Loglikelihood Calculation
def Loglikelihood(z1,z2):
    z10, z20 = 0.954166667,	1.079328859
    # z10, z20 = 0.809361060171, 2.03669874783
    # case 3
    # z10, z20 = 0.809358651978, 2.03723510311
    # z10, z20 = 0.81, 2
    V1 = sig1p*sig1m
    V1e = sig1p - sig1m
    V2 = sig2p * sig2m
    V2e = sig2p - sig2m
    V1f = V1 + V1e*(z1-z10)
    V2f = V2 + V2e * (z2 - z20)
    L2t = 1/(1-p**2)*((z1-z10)**2/V1f-2*p*(z1-z10)*(z2-z20)/np.sqrt(V1f*V2f)+(z2-z20)**2/V2f)
    return L2t

# Grid calculations
xmax, ymax = 4, 2.2
(a, b) = (0.01, 0.01)
x2 = np.arange(0, xmax, a, dtype=np.float)
y2 = np.arange(0.5, ymax, b, dtype=np.float)
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

# read data for official 68% and 95% CL contours
expdata = np.genfromtxt('HIGG-2018-51-Llh-2d-GridFull.txt')
xExp = expdata[:,0]
yExp = expdata[:,1]

# Make validation plot

matplotlib.rcParams['xtick.major.pad'] = 15
matplotlib.rcParams['ytick.major.pad'] = 15

plt.contourf(X2,Y2,Z,[10**(-5),2.3,5.99],colors=['0.5','0.75'])
plt.plot(xExp,yExp,'.',markersize=4, color = 'blue', label="ATLAS official")
plt.plot([1],[1], '*', c='k', ms=6, label="SM")

plt.xlim((0,2))
plt.ylim((0,3))
plt.minorticks_on()
plt.tick_params(labelsize=20, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

# Plot
plt.legend(loc='upper right', fontsize=12)
plt.xlabel(r'$\mu_{bb}^{WH}$', fontsize=20)
plt.ylabel(r'$\mu_{bb}^{ZH}$', fontsize=20)
# plt.title("$\mu$ from ATLAS-HIGG-2016-22 (VGaussian, 68% CL, fitted corr)")
plt.title("$\mu$ from ATLAS-HIGG-2018-51 (VGaussian, 95% CL, fitted corr)")
# plt.title("$\mu$ from ATLAS-HIGG-2016-22 (VGaussian using VGaussian data, corr = -0.27)")
# plt.title("$\mu$ from ATLAS-HIGG-2016-22 (VGaussian, provided 1D, corr = -0.27)")

fig.set_tight_layout(True)

# fig.savefig('mu_ggH-VBF_2D_VGaussian-fitted-corr-68CL.pdf')
# fig.savefig('HIGG-2016-22-mu-VGaussian-fit2d-VGaussian-Fig8a-68CL.pdf')
# fig.savefig('HIGG-2016-21-mu-Poisson-fit1d-VGaussian-AuxFig23.pdf')
# fig.savefig('HIGG-2016-21-mu-VGaussian-fit1d-VGaussian-AuxFig23.pdf')
fig.savefig('ATLAS-HIGG-2018-51 (VGaussian, 95% CL, fitted corr).pdf')

plt.show()