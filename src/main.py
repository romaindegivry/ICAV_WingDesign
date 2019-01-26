import sys, os
print(sys.version)
print(os.getcwd())

#Import modules
import logging as log

import numpy as np
import aero as ae #aerodynamic models
#from . import structural as stc #structural models
#from . import performance as prf #aircraft performance models
import modelling as mdl #model interfaces
#from . import report as rep #reporting

#Constants
global GSI
GSI = 9.81

airfoil = "../in/6043.dat"
wingOptions = {"wRChord" : 0.3,
               "wTChord" : 0.15,
               "span" : 2.,
               "wAngle" : 0.,
               "wAirfoil" : airfoil}

htOptions = {"htRChord" : 0.2,
               "htTChord" : 0.1,
               "htSpan" : 0.7,
               "htAngle" : 0.,
               "htAirfoil" : airfoil,
               "htRLoc" : np.array([[-1.],[0.],[0.]])}

vtOptions = {"vtRChord" : 0.2,
             "vtTChord" : 0.1,
             "vtHeight" : 0.40,
             "vtAirfoil" : airfoil,
             "vtRLoc" : np.array([[-1.],[0.],[0.]])}

meshOptions = {"n_chordwise" : 15,
               "n_spanwise" : 15}

model = mdl.SimpleUAV()
model.aerodynamics(**wingOptions,**htOptions,**vtOptions,**meshOptions)
model._aero.add_case(name='cruise',alpha=+2)
model._aero.run()
model._aero._session.show_geometry()
#begin by importing the airfoil
