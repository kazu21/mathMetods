from matplotlib.pyplot import savefig, plot, show

from numpy import *

"""
===== ==========================================================
Name  Description
===== ==========================================================
Nx    The total number of mesh cells; mesh points are numbered
      from 0 to Nx.
F     The dimensionless number a*dt/dx**2, which implicitly
      specifies the time step.
T     The stop time for the simulation.
I     Initial condition (Python function of x).
a     Variable coefficient (constant).
L     Length of the domain ([0,L]).
x     Mesh points in space.
t     Mesh points in time.
n     Index counter in time.
u     Unknown at current/new time level.
u_1   u at the previous time level.
dx    Constant mesh spacing in x.
dt    Constant mesh spacing in t.
===== ==========================================================
"""


def solver_FE_simple():
    """
    Simplest expression of the computational algorithm
    using the Forward Euler method and explicit Python loops.
    For this method F <= 0.5 for stability.
    """

    L = 1.
    a = 1
    T = 0.1
    F = 0.5
    Nx = 50

    x = linspace(0, L, Nx + 1)  # mesh points in space
    dx = x[1] - x[0]
    dt = F * dx ** 2 / a
    Nt = int(round(T / float(dt)))
    u = zeros(Nx + 1)
    u_1 = zeros(Nx + 1)

    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx + 1):
        u_1[i] = exp(-0.5 * ((x[i] - L / 2.0) ** 2) / 0.05 ** 2)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_1[i] + F * (u_1[i - 1] - 2 * u_1[i] + u_1[i + 1])

        # Insert boundary conditions
        u[0] = 0
        u[Nx] = 0

        # Switch variables before next step
        u_1, u = u, u_1
        if n % 10 == 0:
            plot(x, u)
    show()


solver_FE_simple()
