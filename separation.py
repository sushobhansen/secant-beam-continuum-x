import numpy as np
import matplotlib.pyplot as plt

t,row,col,dt1,dn1,dt2,dn2 = np.loadtxt('separation.txt',unpack=True)

nsteps = t[len(t)-1]
nsteps = int(nsteps)


nx = len(t)/nsteps
nx = int(nx)+1 #No of points
T = np.arange(1,nsteps+1)/nsteps
dt_max = 0.112443640581129*np.ones(nsteps)
dn_max = 0.112443640581129*np.ones(nsteps)
dn = np.zeros((nsteps,nx))


for i in range(nsteps):
	for j in range(nx):
		if j==(nx-1):
			dn[i,j] = dn2[(nx-1)*i+(nx-2)]
		else:
			dn[i,j] = dn1[(nx-1)*i+j]

plt.plot(T,dn[:,nx-1],label='x = 150 mm')
plt.plot(T,dn_max,label='Max Separation')
plt.xlabel('Time')
plt.ylabel('Normal Separation (mm)')
plt.legend()

plt.show()