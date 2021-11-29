import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from n_pnd import *

N = [5, 5] # number of linked pendulums
SLOW_MOTION = 0.25 # slow motion factor in animation
dt = 0.001 # infinitesimal time increment in seconds
T = 100 # length of time to be simulated
G = 9.8
coef = setup_coef(N)

# initial conditions
theta0 = []
omega0 = []
for i in range(len(N)):
    theta0.append(np.tile(2.0, N[i])) # counterclockwise angle relative to the vertical for each rod
    omega0.append(np.tile(0.0, N[i])) # counterclockwise angular velocity for each rod
theta0[1][0] += 1.0e-10

print("INITIAL THETA")
print(theta0[0])
print(theta0[1])

coef = setup_coef(N)

fig, ax = plt.subplots()
time_text = ax.text(-0.9*max(N), 0.9*max(N), "t=0")
data = []
init_data = [np.copy(theta0), np.copy(omega0)]
for i in range(len(N)):
    data.append(evolve_pnd(init_data[0][i], init_data[1][i], T, dt, G, coef))
lines = [ax.plot(data[i][0][0], data[i][1][0])[0] for i in range(len(data))]
draw = lines + [time_text]
nframes = int(len(data[0][0])*dt*SLOW_MOTION/0.02-1)
def animate(i):
    for j in range(len(lines)):
        lines[j].set_xdata(data[j][0][int(i*0.02/dt/SLOW_MOTION)])
        lines[j].set_ydata(data[j][1][int(i*0.02/dt/SLOW_MOTION)])
    time_text.set_text("t="+str(int(i*20)/SLOW_MOTION/1000.0))
    return draw
ani = animation.FuncAnimation(fig, animate, frames=nframes, interval=20, save_count=nframes)
plt.xlim([-1.1*max(N), 1.1*max(N)])
plt.ylim([-1.1*max(N), 1.1*max(N)])

if len(N) >= 2:
    fig2, ax2 = plt.subplots()
    print("Calculating Lyapunov")
    print(theta0)
    data_t_0 = [data[0][2], data[0][3]]
    data_t_1 = [data[1][2], data[1][3]]
    lyap = lyapunov(data_t_0, data_t_1, dt)
    print(lyap[1].shape)
    print(lyap[2])
    ax2.plot(lyap[0], lyap[1])

plt.show()
