import numpy as np

# numerical model parameters

# model size
from matplotlib.pyplot import plot, show, pcolor

xSize = 1e6
ySize = 1.5e6

# number of nodes
xnum = 50
ynum = 30

# grid steps
xstep = xSize / xnum
ystep = ySize / ynum

tnum = 20

# material properties
k = 3  # thermal conductivity
cp = 1e3  # Heat capacity
rho = 3200  # density
rhocp = rho * cp  # volumetric heat capacity

kappa = k / rhocp  # thermal diffusivity

# basic nodes
x = np.arange(0, xSize, xstep)
y = np.arange(0, ySize, ystep)

# timestep
dtExp = min(xstep, ystep) ** 2 / 3 / kappa
dt = 1.5 * dtExp

# temperature profile

tback = 1e3  # background temperature
tWave = 1300  # temperature wave

# create temperature array for implicit solving
t0imp = tback * np.ones((ynum, xnum))
t1imp = t0imp[:]

for i in range(ynum):
    for j in range(xnum):
        if ySize * 0.3 < y[i] < ySize * 0.7 and xSize * 0.3 < x[j] < xSize * 0.7:
            t0imp[i, j] = tWave

# time cycle
for t in range(tnum):
    L = np.zeros((xnum * ynum, xnum * ynum))  # matrics coef init
    R = np.zeros((xnum * ynum, 1))

    for i in range(ynum):
        for j in range(xnum):
            k = j * ynum + i  # global index
            if i == 0 or i == ynum-1 or j == 0 or j == xnum-1:
                if i == 0:
                    L[k, k] = 1
                    R[k, 0] = tback
                if i == ynum - 1:
                    L[k, k] = 1
                    R[k, 0] = tback
                if j == 0 and 0 < i < ynum-1:
                    L[k, k] = 1
                    R[k, 0] = tback
                if j == xnum - 1 and 0 < i < ynum-1:
                    L[k, k] = 1
                    R[k, 0] = tback
            else:
                L[k, k - ynum] = -kappa / xstep ** 2
                L[k, k + ynum] = -kappa / xstep ** 2
                L[k, k - 1] = -kappa / ystep ** 2
                L[k, k + 1] = -kappa / ystep ** 2
                L[k, k] = 1 / dt + 2 * kappa / xstep ** 2 + 2 * kappa / ystep ** 2
                R[k, 0] = t0imp[i, j] / dt
    S = np.linalg.solve(L, R)
    for i in range(ynum):
        for j in range(xnum):
            k = j * ynum + i
            t1imp[i, j] = S[k]

    pcolor(x / 1000, y / 1000, t1imp)

    t0imp = t1imp[:]

show()
