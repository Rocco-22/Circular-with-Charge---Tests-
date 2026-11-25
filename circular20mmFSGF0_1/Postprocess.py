import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def ReadPhi(name):
    L = open(name).readlines().splitlines()
    nx0,Nx,ny0,Ny = [int(x) for x in L[0].split()]
    X = np.array([float(x) for x in L[2:2+Nx].split()]) 
    Y = np.array([float(x) for x in L[2+Nx+1:2+Nx+Ny].split()])
    Phi = np.array([float(x) for x in L[2+Nx+Ny+6:2+Nx+Ny+6+Nx*Ny].split()]).reshape(Ny,Nx)
    return X,Y,Phi

def ReadMaterial(name):
    L = open(name).readlines().splitlines()
    nx,ny,z = [int(x) for x in L[0].split()]



Ex, Ey = -dx, -dy                             # E = -grad(Phi)
X, Y = 1000*X, 1000*Y                         # convert to mm
# Deep-metal nodes are exactly 0 in both components: hide them
mask = (dx == 0) & (dy == 0)
Phi = np.ma.array(Phi, mask=mask) 
Ex = np.ma.array(Ex, mask=mask)
Ey = np.ma.array(Ey, mask=mask)