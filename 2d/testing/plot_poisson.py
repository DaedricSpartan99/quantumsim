import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

"""
    File structure

    x y numerical exact difference
"""
def readSolution():
    U = []
    x = []
    y = []
    with open("poisson.dat") as f:
        lines = f.readlines()

        for l in lines:
            line = l.split(' ')
            x.append(float(line[0]))
            y.append(float(line[1]))
            U.append(float(line[2]))

        return x, y, U

"""
    Convergence data structure

    N h |U - U_h|
"""
def readConvergence():
    N = []
    h = []
    diff_norm = []
    with open("poisson_convergence.dat") as f:
        lines = f.readlines()
        
        for l in lines:
            line = l.split(' ')
            N.append(int(line[0]))
            h.append(float(line[1]))
            diff_norm.append(float(line[2]))

        return N, h, diff_norm


# plot numerical solution
X, Y, U = readSolution()

X = np.array(X)
Y = np.array(Y)
U = np.array(U)

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
surf = ax.plot_trisurf(X, Y, U, cmap=cm.jet, linewidth=0.1)
fig.colorbar(surf, shrink=0.5, aspect=5)

ax.set_title('Numerical solution')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('U')


# plot convergence graph
N, h, diff = readConvergence() 

h = np.array(h)
diff = np.array(diff)

p = np.polyfit(h, diff, 1)
fit = lambda x: p[0] * x + p[1]

plt.figure()
plt.scatter(h, diff, label="Data")
plt.plot(h, fit(h), label="Linear fit")

plt.grid()
plt.legend()

plt.title("Convergence test")
plt.xlabel("h")
plt.ylabel(r"|U - U_h|")

plt.show()
