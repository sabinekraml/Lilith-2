#import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

## only for example, use your grid
#z = np.linspace(0, 1, 15)
#x = np.linspace(0, 1, 15)
#y = np.linspace(0, 1, 15)

#X, Y, Z = np.meshgrid(x, y, z)

## Your 4dimension, only for example use yours
#U = np.exp(-(X/2) ** 2 - (Y/3) ** 2 - Z ** 2)

## Creating figure
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

## Creating plot
#ax.scatter3D(X, Y, Z, c=U, alpha=0.7, marker='.')
#plt.show()

#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#import numpy as np

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

#x = np.random.standard_normal(100)
#y = np.random.standard_normal(100)
#z = np.random.standard_normal(100)
#c = np.random.standard_normal(100)

#print("c = ",c)

#img = ax.scatter(x, y, z, c=c, cmap=plt.hot())
#fig.colorbar(img)
#plt.show()

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# domains

data = np.genfromtxt("/home/Willy/Lilith/Lilith-2/validations/STU/STU.out")
x = data[:,0]
y = data[:,1]
z = data[:,2]
print("test")
#x = np.logspace(-1.,np.log10(5),50) # [0.1, 5]
#y = np.linspace(6,9,50)             # [6, 9]
#z = np.linspace(-1,1,50)            # [-1, 1]

# convert to 2d matrices
Z = np.outer(z.T, z)        # 50x50
X, Y = np.meshgrid(x, y)    # 50x50

# fourth dimention - colormap
# create colormap according to x-value (can use any 50x50 array)
color_dimension = X # change to desired fourth dimension
minn, maxx = color_dimension.min(), color_dimension.max()
norm = matplotlib.colors.Normalize(minn, maxx)
m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
m.set_array([])
fcolors = m.to_rgba(color_dimension)

# plot
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X,Y,Z, rstride=1, cstride=1, facecolors=fcolors, vmin=minn, vmax=maxx, shade=False)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
