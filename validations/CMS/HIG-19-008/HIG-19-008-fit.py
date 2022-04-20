from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize as optimize
import sys

# Open the 1D grid files

# Choose VBF - ggH - ttH - VH
f = open('validations/CMS/HIG-19-008/HIG-19-008-muttHmutH-Grid95.txt', 'r')

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
    z10, z20 = 5.7, 0.92

    V1 = sig1p * sig1m
    V1e = sig1p - sig1m
    V2 = sig2p * sig2m
    V2e = sig2p - sig2m
    V1f = V1 + V1e * (z1 - z10)
    V2f = V2 + V2e * (z2 - z20)
    L2t = 1 / (1 - p ** 2) * ( (z1 - z10) ** 2 / V1f - 2 * p * (z1 - z10) * (z2 - z20) / np.sqrt(V1f * V2f) + (z2 - z20) ** 2 / V2f )
#    print("L2t = ", L2t)
    return L2t

def fit5para(xr, yr):
    A = []
    B = []
    for i in range(0, len(xr)):
        A.append([xr[i], yr[i]])
        B.append(5.99)
    A = np.array(A)
    B = np.array(B)
    guess = (4.1, 4, 0.21, 0.2, 0.2)
#    guess = (0.6, 0.5, 0.6, 0.5, 0.3)
    AR, pcov = optimize.curve_fit(func, A, B, guess)
    return AR


# Find parameters
def func7(X, mu1, mu2, sig1p, sig1m, sig2p, sig2m, p):
    z1, z2 = X[:,0], X[:,1]
    z10, z20 = mu1, mu2

    V1 = sig1p * sig1m
    V1e = sig1p - sig1m
    V2 = sig2p * sig2m
    V2e = sig2p - sig2m
    V1f = V1 + V1e * (z1 - z10)
    V2f = V2 + V2e * (z2 - z20)
    L2t = 1 / (1 - p ** 2) * ( (z1 - z10) ** 2 / V1f - 2 * p * (z1 - z10) * (z2 - z20) / np.sqrt(V1f * V2f) + (z2 - z20) ** 2 / V2f )
#    print("L2t = ", L2t)
    return L2t

def fit7para(xr, yr):
    A = []
    B = []
    for i in range(0, len(xr)):
        A.append([xr[i], yr[i]])
        B.append(5.99)
    A = np.array(A)
    B = np.array(B)
    guess = (5.7, 0.92, 4.1, 4, 0.21, 0.2, 0.2)
    AR, pcov = optimize.curve_fit(func7, A, B, guess)
    return AR


#ff = fit5para(x, y)
#sig1p, sig1m, sig2p, sig2m, p = ff[0], ff[1], ff[2], ff[3], ff[4]

#print("\np =", p)
#print("\nParameters:")
#print("sig1p, sig1m, sig2p, sig2m, corr =", ff[0], ",", ff[1], ",", ff[2], ",", ff[3], ",", ff[4])


ff = fit7para(x, y)
mu1, mu2, sig1p, sig1m, sig2p, sig2m, p = ff[0], ff[1], ff[2], ff[3], ff[4], ff[5], ff[6]

print("\np =", p)
print("\nParameters:")
print("mu1, mu2, sig1p, sig1m, sig2p, sig2m, corr =", ff[0], ",", ff[1], ",", ff[2], ",", ff[3], ",", ff[4], ",", ff[5], ",", ff[6])


#sig1p, sig1m, sig2p, sig2m, corr = 4.117334822885522 , 4.012796468185626 , 0.2618919524977251 , 0.2368266075188232 , -0.25111013734589455
#mu1, mu2, sig1p, sig1m, sig2p, sig2m, corr = 5.983051099311637 , 0.9279559177987796 , 4.0512447344911395 , 4.0757276694986855 , 0.2598745711648347 , 0.23901752974630947 , -0.2508891023877273
