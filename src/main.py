import sys, os
print(sys.version)
print(os.getcwd())

#Import modules
import logging as log

import numpy as np
import aero as ae #aerodynamic models
#from . import structural as stc #structural models
#from . import performance as prf #aircraft performance models
#from . import modelling as mdl #model interfaces
#from . import report as rep #reporting

#Constants
global GSI
GSI = 9.81
aero = ae.Aero()
wing = ae.Wing()
wing.autoSections(0.3,0.15,1.,"../in/6043.dat")
wing.autoWing(15, 15)
aero.set_wing(wing)
aero.build_geometry(1.,1.,1.)
aero.add_case(name='cruise',alpha=+2)
aero.run()
print(aero.results)
#begin by importing the airfoil
