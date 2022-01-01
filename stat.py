import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

CSV = "dat.csv"
RMV = 0

def list_float(l): # converts a list of some type to a list of floats
    ret = []
    for i in range(len(l)):
        ret.append(float(l[i]))
    return ret
def linear_func(ydata): # function intended to linearize data
    return np.exp(np.asarray(ydata))
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
    print(pearson[i])

colors = ['C0', 'C1', 'C2', 'y']
plt.figure(0)
for i in range(len(data)):
    plt.plot(list_float(N), data[i], label=angles[i], color=colors[i%len(colors)])
plt.xlabel("Number of pendulums")
plt.ylabel("Lyapunov exponent $\\lambda$ (s${}^{-1}$)")
plt.legend(title="Angle (rad)")

plt.figure(1)
for i in range(len(data)):
    plt.plot(list_float(N), linear_func(data[i]), label=angles[i], color=colors[i%len(colors)])
    plt.plot(list_float(N), regress_func(pearson[i]), color=colors[i%len(colors)], linestyle='--')
plt.title("Linear regression model for increasing chaos in connected pendulums")
plt.xlabel("Number of pendulums")
plt.ylabel("exp$( \\frac{1}{\\lambda} )$")

fig, ax = plt.subplots(2, 2)
fig.suptitle("Standardized residual plots")
residuals = []
fitted = []
for i in range(len(data)):
    residuals.append(linear_func(data[i])-regress_func(pearson[i]))
    fitted.append(linear_func(data[i]))
    ax[int(i/2),i%2].scatter(fitted[i], residuals[i]/np.std(residuals[i]), color=colors[i%len(colors)])

fig2, ax2 = plt.subplots(2, 2)
for i in range(len(residuals)):
    stats.probplot(residuals[i], plot=ax2[int(i/2),i%2])
    ax2[int(i/2),i%2].set_title("")

plt.show()

