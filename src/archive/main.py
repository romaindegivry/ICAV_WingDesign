import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

import deps.beamfe.beamfe as FEM
GSI = 9.81

#problem properties
fibre = {'E' : 33e9, 'rho' : 1666.6,'sigTens': 400,'sigComp' : 200}
foam = {'E' : 1.9e9, 'rho' : 30.0,'sigTens': 22,'sigComp' : 22}

def stressArray(rot,E,y,dx):
    #differentiate the rotation array
    d2w = np.diff(rot)/dx
    #compute the stress and return
    try:
        return -E*y[1:]*d2w
    except IndexError:
        return -E*y*d2w



class Material:

    def __init__(self,E,rho):
        self.E = E
        self.rho = rho

    def set_failure(self,sigComp,sigTens):
        self.sigTens = sigTens
        self.sigComp = sigComp

class Wing:
    """
    This class represents a wing design whith all associated properties
    and available calculations
    """

    def __init__(self,xStart,xEnd,Nnodes = 100):
         self.L = np.linalg.norm(xStart-xEnd) #Length of the beam
         self.normal = (xEnd-xStart)/self.L
         self.start = xStart
         self.Nnodes = Nnodes
         self.nodes = np.linspace(0,self.L,num=Nnodes)

    def set_material(self, fibre, foam):
        self.fibre = fibre
        self.foam = foam

    def EI(self,b,t):
        #use formula to compyte the equivalent EI
        Ef = self.fibre.E
        Ec = self.foam.E
        c = self.c
        EIz = Ef*((b*t**3)/(6) + (b*t*(self.c+t)**2)/2) + Ec*(b*self.c**3)/12
        EIy = Ef*((b**3*t)/6) + Ec*(self.c*b**3)/12
        return EIy, EIz

    def EA(self,b,t):
        #compute the area of the plate
        c = self.c
        return self.fibre.E *2*t*b + self.foam.E *c*b

    def density(self,b,t):
        c = self.c
        return (2*t*self.fibre.rho + c*self.fibre.rho)/(c+2*t)

    def set_chord_distrib(self,c):
        self.c = c

    def buildLoadCase(self,m,b,t,safetyFactor=1.2,n=2.5,LD_2D=15):
        c = self.c
        #Load cases
        weight_N = m*GSI*safetyFactor*n
        self.load = np.zeros((6,Nnodes))
        lift = (weight_N*2)/(np.pi*self.L**2)*np.sqrt(self.L**2-self.nodes**2)
        drag = lift/LD_2D
        weight = -GSI*self.density(b,t)*b*(2*t+c)*n*safetyFactor
        self.load[2,:] = lift + weight
        self.load[1,:] = drag
        assert(np.abs(np.trapz(lift,x=self.nodes)/(GSI*2.5*1.2)*2 - 6.9) < 0.05) #check that lift is below 50g off

    def solve(self,b,t):
        load = FEM.interleave(self.load.T)
        EIy,EIz = self.EI(b,t)
        prob = FEM.BeamFE(self.nodes,self.density(b,t),self.EA(b,t),EIy,EIz)
        prob.set_boundary_conditions(left='C',right='F')
        prob.set_dofs([True,True,True,False,True,True])
        return prob.static_deflection(prob.distribute_load(load))[0]

    def failure(self,b,t):
        c = self.c
        x = self.solve(b,t)
        δ = FEM.deleave(x)
        disp = {'x' :  δ[0,:],
                'y' :  δ[1,:],
                'z' :  δ[2,:],
                'rx' : δ[3,:],
                'ry' : δ[4,:],
                'rz' : δ[5,:]}

        #failure ratio
        fibRatio = (np.abs(stressArray(disp['rz'],self.fibre.E,c/2+t,self.nodes[1]-self.nodes[0]))
                  + np.abs(stressArray(disp['ry'],self.fibre.E,b/2  ,self.nodes[1]-self.nodes[0])))/min(self.fibre.sigTens,self.fibre.sigComp)

        fomRatio = (np.abs(stressArray(disp['rz'],self.foam.E,c/2,self.nodes[1]-self.nodes[0]))
                  + np.abs(stressArray(disp['ry'],self.foam.E,b/2,self.nodes[1]-self.nodes[0])))/min(self.foam.sigTens,self.foam.sigComp)
        return max(np.max(fomRatio),np.max(fibRatio))

    def mass(self,b,t):

        return np.trapz(self.density(b,t)*b*(self.c+2*t),x=self.nodes)


#define materials
fibre = Material(33e9,1666.6); fibre.set_failure(400e6,200e6)
foam = Material(1.9e9,30.0); foam.set_failure(400e6,200e6)
Nnodes =100
root = np.array([0,0,0])
tip = np.array([0,1996.0/2*1e-3,0])

#initialise the problem
halfWing = Wing(root,tip, Nnodes)
halfWing.set_material(fibre,foam)
halfWing.set_chord_distrib(np.linspace(0.03,0.015,num=Nnodes))


#define optimisation problem
def optObj(x,*args):
    print("Function evaluation")
    b = x[0]
    t = x[1:Nnodes+1]
    halfWing.buildLoadCase(6.9,b,t)
    #add linear penalty
    fail = halfWing.failure(b,t)
    if fail > 1:
        penalty = fail
        print("fail")
    else:
        penalty = 0
    print("objective", halfWing.mass(b,t) + penalty*1e6)
    return halfWing.mass(b,t) + penalty*1e6

#run optimisation code
#initial conditions
b0 = 1e-2
t0 = np.ones(Nnodes)*1e-2
x0 = np.hstack((b0,t0))
bounds = [(x0[i]*0.1,x0[i]*10) for i in range(len(x0))]
print("Begin optimisation")
res = scipy.optimize.differential_evolution(optObj,bounds,maxiter=10)

print(res)
