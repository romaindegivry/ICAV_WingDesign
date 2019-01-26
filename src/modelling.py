import numpy as np

import aero as ae

"""
This sub package is in charge of keeping the connections between subsystem
models of the aircraft. The Model class contains baseline functionality and
will then be extented to build specific models
"""


class Model:
    """Contains bindings between models"""

    def __init__(self):
        self.aeroMdl = None
        self.structMdl = None
        self.components = {}
        self.powerMdl = None

    def add_Aero(aeroMdl):
        #a bit of protection
        if not isinstance(aeroMdl, ae.Aero):
            raise ValueError("Aerodynamic model must be instance of ae.Aero")
        self.aeroMdl = aeroMdl

class SimpleUAV(Model):
    """Class that builds baselime models for a simple (HT,VT, Wing) uav"""

    def __init__(self, definition= None):
        """Creates the class but may require"""
        if definition == None:
            return
        else:
            #actually initialise
            #Here we raise an exception until this is implemented
            raise NotImplementedError()

    def aerodynamics(self, wRChord = 0.3,
                           wTChord = 0.15,
                           span = 2.,
                           wAngle = 0.0,
                           wAirfoil = "",
                           htRChord= 0.2,
                           htTChord= 0.1,
                           htSpan = 1.,
                           htAngle = 0.0,
                           htAirfoil = "",
                           htRLoc = np.array([[-1.],[0.],[0.]]),
                           vtRChord = 0.2,
                           vtTChord = 0.1,
                           vtHeight = 0.5,
                           vtAirfoil = "",
                           vtRLoc = np.array([[-1.],[0.],[0.]]),
                           n_chordwise = 15,
                           n_spanwise = 15,
                           offset=np.zeros((3,1))):
        """Creates an aircraft model at offset"""
        #mean aerodynamic chord
        self.MAC = (wRChord + wTChord)/2.
        self.c = self.MAC

        #Reference span
        self.bref = span

        #reference area
        self.Aref = span*self.MAC
        #Wing
        self._wing = ae.Wing()
        self._wing.autoSections(wRChord,wTChord,span/2,
                                wAirfoil,rootPos =offset)
        self._wing.autoWing(n_chordwise,n_spanwise)

        #Horizontal tail
        self._htail = ae.HTail()
        self._htail.autoSections(htRChord, htTChord, htSpan/2, htAirfoil,\
                            rootPos = htRLoc+offset)
        self._htail.autoHTail(n_chordwise,n_spanwise)

        #vertical tail
        self._vtail = ae.VTail()
        self._vtail.autoSections(vtRChord, vtTChord, vtHeight, vtAirfoil,\
                            rootPos = vtRLoc+offset)
        self._vtail.autoVTail(n_chordwise,n_spanwise)

        #build the model
        self._aero = ae.Aero()
        self._aero.set_wing(self._wing)
        self._aero.set_htail(self._htail)
        self._aero.set_vtail(self._vtail)
        self._aero.build_geometry(self.Aref, self.MAC, self.bref)



if __name__ == "__main__":
    #if the submodule is executed perform tests
    import unittest

    class testAeroModel(unittest.TestCase):

        def test_init(self):
            model = SimpleUAV()

        def test_aerodynamics(self):
            model = SimpleUAV()
            model.aerodynamics()
    unittest.main(verbosity=2)
