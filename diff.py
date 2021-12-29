import numpy as np
import matplotlib.pyplot as plt
from n_pnd import *

n_rods = 200
angles = [2.0]
N = []
for i in range(len(angles)):
    N.append(n_rods)
    N.append(n_rods)
dt = 0.001 # infinitesimal time increment in seconds
T = 20 # length of time to be simulated
L = 50.0
G = 9.8/L # effective gravity
coef = setup_coef(N)

print("Testing n="+str(n_rods))

# initial conditions
theta0 = []
omega0 = []
for i in range(len(angles)*2):
    theta0.append(np.tile(angles[int(i/2.0)], N[i])) # counterclockwise angle relative to the vertical for each rod
    omega0.append(np.tile(0.0, N[i])) # counterclockwise angular velocity for each rod
    if i % 2 == 1:
        theta0[i][0] += 1.0e-15

coef = setup_coef(N)

data = []
init_data = [np.copy(theta0), np.copy(omega0)]
for i in range(len(N)):
    data.append(evolve_pnd(init_data[0][i], init_data[1][i], T, dt, G, coef))
    print("Simulated " + str(i+1)+"/"+str(len(N)))
print("")

x0 = []
y0 = []
x1 = []
y1 = []
for i in range(int(len(N)/2)):
    x0.append(data[i*2][0])
    y0.append(data[i*2][1])
    x1.append(data[i*2+1][0])
    y1.append(data[i*2+1][1])
x0_0 = np.asarray(x0)[:,:,1]
y0_0 = np.asarray(y0)[:,:,1]
x1_0 = np.asarray(x1)[:,:,1]
y1_0 = np.asarray(y1)[:,:,1]
d_0 = np.sqrt(np.square(x1_0-x0_0) + np.square(y1_0-y0_0)).T
plt.plot(np.arange(dt, T, dt), np.log(d_0[1:]))
x0_5 = np.asarray(x0)[:,:,int(n_rods/2)]
y0_5 = np.asarray(y0)[:,:,int(n_rods/2)]
x1_5 = np.asarray(x1)[:,:,int(n_rods/2)]
y1_5 = np.asarray(y1)[:,:,int(n_rods/2)]
d_5 = np.sqrt(np.square(x1_5-x0_5) + np.square(y1_5-y0_5)).T
plt.plot(np.arange(dt, T, dt), np.log(d_5[1:]))
x0_1 = np.asarray(x0)[:,:,n_rods]
y0_1 = np.asarray(y0)[:,:,n_rods]
x1_1 = np.asarray(x1)[:,:,n_rods]
y1_1 = np.asarray(y1)[:,:,n_rods]
d_1 = np.sqrt(np.square(x1_1-x0_1) + np.square(y1_1-y0_1)).T
plt.plot(np.arange(dt, T, dt), np.log(d_1[1:]))

#plt.savefig("bah/diff"+str(n_rods)+".png")
plt.show()
