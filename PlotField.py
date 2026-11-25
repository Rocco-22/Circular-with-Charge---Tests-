import numpy as np
import matplotlib.pyplot as plt

def ReadFile(name):
    d = np.loadtxt(name)
    return d[0, 1:], d[1:, 0], d[1:, 1:]     # X, Y, M[j, i]

X, Y, Phi = ReadFile('PhiMatrix.dat')
_, _, dx  = ReadFile('dxPhiMatrix.dat')
_, _, dy  = ReadFile('dyPhiMatrix.dat')

Ex, Ey = -dx, -dy                             # E = -grad(Phi)
X, Y = 1000*X, 1000*Y                         # convert to mm
# Deep-metal nodes are exactly 0 in both components: hide them
mask = (dx == 0) & (dy == 0)
Phi = np.ma.array(Phi, mask=mask) 
Ex = np.ma.array(Ex, mask=mask)
Ey = np.ma.array(Ey, mask=mask)

fig, ax = plt.subplots(figsize=(7, 7))
ax.contour(X, Y, Phi, 20, colors='black', linewidths=0.5)
ax.streamplot(X, Y, Ex, Ey, density=0.6, color='b', linewidth=0.7,broken_streamlines=False,arrowstyle='-')
ax.set_xlabel('x [mm]'); ax.set_ylabel('y [mm]')
ax.ticklabel_format(style='plain', axis='both', scilimits=(0, 0))
plt.savefig('fieldlines.png', dpi=200)