from math import cos as cos
from math import pi as pi
from matplotlib import pyplot as plt
import numpy as npy

#regional isostazy model
l = 1e6
hm = 5e3
roc = 26e2
rom = 32e2
g = 9.81
lmb = l / 2
d = 1e24

n = 1001
dx = l / (n - 1)

x = [i for i in range(int(-l / 2), int(l / 2) + 1, int(dx))]

h = [hm * (1 + cos(2 * pi * x[i] / lmb)) / 2 if x[i] >= -lmb / 2 and x[i] <= lmb / 2 else 0 for i in range(n)]

pc = roc * g * (dx ** 4) / d
pm = rom * g * (dx ** 4) / d

a = npy.matrix(npy.empty(shape=(n,n)))
b = npy.array(npy.empty(n))

a[0,0] = 1
a[1,1] = 1
a[n - 2,n - 2] = 1
a[n - 1,n - 1] = 1
b[0] = 0
b[1] = 0
b[n - 2] = 0
b[n - 1] = 0

for i in range(2,n - 2):
    a[i,i] = 6 + pm
    a[i,i + 1] = -4
    a[i,i - 1] = -4
    a[i,i + 2] = 1
    a[i,i - 2] = 1
    b[i] = pc * h[i]

w = npy.linalg.solve(a,b)

r = [h[i] - w[i] for i in range(n)]

plt.title("region")
plt.plot(range(n), [-i for i in w], 'r--', label="w")
plt.plot(range(n), r, 'b', label="r")
plt.plot(range(n), h, 'g', label="h")
plt.legend()
plt.show()
