import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

CSV = "dat10.csv"
RMV = 0

def list_float(l):
    ret = []
    for i in range(len(l)):
        ret.append(float(l[i]))
    return ret
def linear_func(ydata): # function intended to linearize data
    return np.exp(1.0/np.asarray(ydata))
def regress_func(p): # regression line from Pearson calculation
    return p.slope*np.asarray(list_float(N))+p.intercept

file = open(CSV, 'r')
temp = file.read().split("\n")
file.close()
data = []
N = temp[0].split(",")
N = N[1:len(N)-RMV]
angles = []
for i in range(1, len(temp)):
    line = temp[i].split(",")
    data.append(list_float(line[1:len(line)-RMV]))
    angles.append(line[0])

spearman = []
pearson = []
for i in range(len(angles)):
    print(angles[i])
    spearman.append(stats.spearmanr(list_float(N), data[i]))
    print(spearman[i])
    pearson.append(stats.linregress(list_float(N), linear_func(data[i])))

colors = ['C0', 'C1', 'C2', 'y']
plt.figure(0)
for i in range(len(data)):
    plt.plot(N, data[i], label=angles[i], color=colors[i%len(colors)])
plt.title("Decreasing chaos for increasing number of connected pendulums")
plt.xlabel("Number of pendulums")
plt.ylabel("Lyapunov exponent $\\lambda$ (s${}^{-1}$)")
plt.legend(title="Angle (rad)")

plt.figure(1)
for i in range(len(data)):
    plt.plot(N, linear_func(data[i]), label=angles[i], color=colors[i%len(colors)])
    plt.plot(N, regress_func(pearson[i]), color=colors[i%len(colors)], linestyle='--')
plt.title("Linear regression model for decreasing chaos in connected pendulums")
plt.xlabel("Number of pendulums")
plt.ylabel("exp$( \\frac{1}{\\lambda} )$")

fig, ax = plt.subplots(2, 2)
fig.suptitle("Standardized residual plots")
for i in range(len(data)):
    residuals = linear_func(data[i])-regress_func(pearson[i])
    ax[int(i/2),i%2].scatter(list_float(N), residuals/np.std(residuals), color=colors[i%len(colors)])

plt.show()

