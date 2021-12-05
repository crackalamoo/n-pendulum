import numpy as np
from n_pnd import *

n_rods = 2
angles = [1.5, 2.0, 2.5, 3.0]
N = []
for i in range(len(angles)):
    N.append(n_rods)
    N.append(n_rods)
SLOW_MOTION = 0.25 # slow motion factor in animation
dt = 0.001 # infinitesimal time increment in seconds
T = 30 # length of time to be simulated
G = 9.8
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

for i in range(len(angles)):
    print("Calculating Lyapunov")
    print(theta0[i*2][0])
    data_t_0 = [data[i*2][2], data[i*2][3]]
    data_t_1 = [data[i*2+1][2], data[i*2+1][3]]
    lyap = lyapunov(data_t_0, data_t_1, dt)
    print(lyap[2])
    print(lyap[3])
    print("")
