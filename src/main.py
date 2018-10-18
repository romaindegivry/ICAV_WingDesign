import numpy as np

import deps.beamfe.beamfe as FEM

GSI = 9.81
#Step 1: Build the discretisation for the beam
def nodes(xStart,xEnd,number):
    #Convert two 3-vectors, xStart, xend into a vector and mesh
    L = np.linalg.norm(xStart-xEnd) #Length of the beam
    n = (xEnd-xStart)/L
    return n, L, np.linspace(0,L,num=number)

Nnodes =50
root = np.array([0,0,0])
tip = np.array([0,1,0])
direction, length, nodes = nodes(root,tip,Nnodes)

#Step 2: Build beam properties
fibre = {'E' : 33e9, 'rho' : 1666.6}
foam = {'E' : 17e9, 'rho' : 30.0}

def EILayered(plate, core, b, t, c):
    #use formula to compyte the equivalent EI
    EIz = plate['E']*((b*t**3)/(6) + (b*t*(c*t)**2)/2) + core['E']*(b*c**3)/12
    EIy = plate['E']*((b**3*t)/6) + core['E']*(c*b**3)/12
    return EIy, EIz

def EALayered(plate, core, b, t, c):
    #compute the area of the plate
    return plate['E']*2*t*b + core['E']*c*b

def densityLayered(plate, core, b,t,c):
    return (2*t*plate['rho'] + c*core['rho'])/(c+2*t)

b = np.ones(Nnodes)*0.025
t = 0.5e-3
c = np.linspace(0.05,0.02,num=Nnodes)

EIy, EIz = EILayered(fibre, foam, b,t,c)
EA = EALayered(fibre, foam, b,t,c)
density = densityLayered(fibre, foam, b,t,c)

#Step 3: Compute the load cases
safetyFactor = 1.2
maxLoadFactor = 2.5
weight_N = 6.9*GSI*safetyFactor

constant_lift = -weight_N/length*np.ones(Nnodes)
load = np.zeros((6,Nnodes))
load[2,:] = constant_lift

#build load vector
load = FEM.interleave(load.T)
#Step 4: Initialise the solver
prob = FEM.BeamFE(nodes,density,EA,EIy,EIz)

prob.set_boundary_conditions('C','F')
prob.set_dofs([True,True,True,False,False,False])
prob.distribute_load(load)
