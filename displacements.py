import numpy as np
import matplotlib.pyplot as plt

xb,yb = np.loadtxt('extmesh.txt',unpack=True)
dispb = np.loadtxt('boundary-displacements.txt')
xi,yi = np.loadtxt('intmesh.txt',unpack=True)
dispi = np.loadtxt('internal-displacements.txt')


dispb = dispb[len(dispb)-2] #One time step earlier than final one to sync with dispi
dispb = dispb[1:]
dispi = dispi[len(dispi)-1]
dispi = dispi[1:]

ub = dispb[0::3]
vb = dispb[1::3]
ui = dispi[0::2]
vi = dispi[1::2]

amp = 10

plt.scatter(xb,yb,color='black',label='Undeformed beam')
plt.scatter(xi,yi,color='red',label='Undeformed element')
plt.scatter(xb+amp*ub,yb+amp*vb,color='orange',marker='+',label='Deformed beam (10x)')
plt.scatter(xi+amp*ui,yi+amp*vi,marker='+',color='blue',label='Deformed element (10x)')
plt.legend()
plt.show()