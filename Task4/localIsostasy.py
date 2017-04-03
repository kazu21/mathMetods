from math import cos as cos
from math import pi as pi
from matplotlib import pyplot as plt

#local isostazy model
l=1e6
hm=5e3
roc=26e2
rom=32e2
g=9.81
lmb=l/2
d=1e24

n=1001
dx=l/(n-1)

x=[i for i in range(int(-l/2), int(l/2)+1, int(dx))]

h=[hm*(1+cos(2*pi*x[i]/lmb))/2 if x[i]>= -lmb/2 and x[i] <= lmb/2 else 0 for i in range(n)]

w=[roc/rom * h[i] for i in range(n)]
r=[h[i]-w[i] for i in range(n)]

plt.title("local")
plt.plot(range(n), [-i for i in w], 'r--', label="w")
plt.plot(range(n), r, 'b', label="r")
plt.plot(range(n), h, 'g', label="h")
plt.legend()
plt.show()
