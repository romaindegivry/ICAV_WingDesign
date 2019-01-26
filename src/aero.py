import numpy as np
from deps.avlwrapper import (Geometry, Surface, Section,  FileAirfoil, Control,
                                       Point, Spacing, Session, Case, Parameter,
                                       ParameterSweep, ProfileDrag)
class Wing:
    """This class contains information about the details of
    the wing aerodynamic model.
    Namely:
     - Geometry information
     - Meshing info
    The root and tip sections of the wing should be specified.
    This class describes multiple ways to build the wing geometry
    Wings are assumed to be symmetric.
    """
    def __init__(self):
        #Declare some variables
        self.rootSection = None
        self.tipSection = None
        self.wingSurface = None

    def autoSections(self, croot, ctip, semispan, airfoilLoc,
                     rootPos = np.zeros((3,1)),sectionKwargs={}):
        """
        rootPos: is the absolute position of the root QUARTER CHORD point
        sectionKwargs: is used to pass keyword arguments to alv,
            check deps.avlwrapper.geometry.Section for details.
            Keyword airfoil will be overwritten.
        """
        #read Airfoil
        AF = FileAirfoil(airfoilLoc)
        airfoil = {"airfoil" : AF}
        #build root
        rootLE = Point(rootPos[0,0]+croot/4.,rootPos[1,0],rootPos[2,0])
        self.rootSection = Section(rootLE,croot,**dict(sectionKwargs,**airfoil))
        #build tip
        tipLE = Point(rootPos[0,0]+ctip/4.,rootPos[1,0]+semispan,rootPos[2,0])
        self.tipSection = Section(tipLE, ctip,**dict(sectionKwargs,**airfoil))


    def autoWing(self,n_spanwise,n_chordwise,name="Wing"):
        """This function builds the wing geometry
        representation from root and tip sections defs"""

        if (self.rootSection == None) or (self.tipSection == None):
            raise NameError("rootSection or tipSection not defined, run autoSections")
        #create the Surface
        self.wingSurface = Surface(name=name,
                                   n_chordwise=n_chordwise,
                                   chord_spacing=Spacing.cosine,
                                   n_spanwise=n_spanwise,
                                   span_spacing=Spacing.cosine,
                                   y_duplicate=0.0, #symmetric
                                   sections=[self.rootSection,self.tipSection])


    def set_root_section(self, croot, airfoilLoc, position):
        """This method overrides autoSections functionality, to let one
        test out small modifications to the assumed wing structure"""
        self.rootSection = None

class AvlSolution:
    """Class to make post processed results accessible, they keep a copy of the
    solution file handy"""

    def __init__(self,resultDict):
        self._raw = resultDict
        #create Output
        self.CL = resultDict
        self.CD = resultDict
        self.

class Aero:
    """This class contains the whole data regarding the aerodyncamic model for
    the aircraft.
    """

    def __init__(self):
        self.surfaces = {"wing" : None,
                         "htail" : None,
                         "vtail" : None}
        self.geometry = None
        self.cases = {}
        self.results = {}

    def set_wing(self,wing):
        """Just sets the wing parameter"""
        self.wing = wing
        self.surfaces["wing"] = wing.wingSurface
        print(wing.wingSurface)

    def build_geometry(self,Aref,cref,bref,name="Aircraft"):
        """Create the geometry for the case"""
        surfs =[]
        #go through surfaces and add them to the geometry
        for i in self.surfaces:
            if self.surfaces[i] != None:
                surfs.append(self.surfaces[i])
        if len(surfs) == 0:
            raise ValueError("Not enough surfaces to build geometry")
        self.geometry = Geometry(name="Aircraft",
                                 reference_area=Aref,
                                 reference_chord=cref,
                                 reference_span=bref,
                                 reference_point=Point(0, 0, 0),
                                 surfaces=surfs)
    def add_case(self,name=None, case=None,**kwargs):
        """Appends a case to the problem, pass a name and kwargs or pre-build case
        through kwarg "case" """
        if case != None and name==None:
            self.cases[case.name] = case
        elif name != None:
            #simple mapping
            self.cases[name] = Case(name = name, **kwargs)
        else:
            raise ValueError("name or case must be supplied")
        pass

    def run(self,preserve_input_file = None : str):
        """Creates the avl session, extracts all info.
        This function requires a list of cases to be specified.
        The results are then saved as a results object for each case
        (allows for easy querying).
        The input file may be kept at the specified address"""

        #extract cases for processing
        allCases = [self.cases[i] for i in self.cases]
        #Check we can actually run this
        if len(allCases) == 0:
            raise ValueError("Not enough cases specified")
        #run, difficut to gurantee avl does not fail
        self._session = Session(geometry=self.geometry, cases=allCases)
        results = self._session.get_results()
        #post process
        for case in results:
            self.results["case"] = AvlSolution(results[case])
        if preserve_input_file != None:
            self._session._write_geometry(name=preserve_input_file)


#add some testing
if __name__ == "__main__":
    #if the submodule is executed perform tests
    import unittest

    class testWingMethods(unittest.TestCase):
        """Test case to check nothing's broken"""
        #share some data
        afLoc = "../in/test_data/test_airfoil.dat"
        rchord = 0.3
        tchord = 0.15
        bspan = 1.

        def test_init(self):
            self.wing = Wing()

        def test_auto_sections(self):
            wing = Wing()
            wing.autoSections(self.rchord, self.tchord, self.bspan, self.afLoc)
            self.assertEqual(wing.rootSection.chord,self.rchord)
            self.assertEqual(wing.tipSection.chord,self.tchord)

        def test_auto_wing(self):
            wing = Wing()
            wing.autoSections(self.rchord, self.tchord, self.bspan, self.afLoc)
            wing.autoWing(15, 15)
            self.assertTrue(isinstance(wing.wingSurface,Surface))

        #set root sections probably does not need a test
    class testAeroMethods(unittest.TestCase):
        afLoc = "../in/test_data/test_airfoil.dat"
        rchord = 0.3
        tchord = 0.15
        bspan = 1.

        def test_init(self):
            Aero()
        def test_setWing(self):
            self.aero = Aero()
            self.wing = Wing()
            self.wing.autoSections(self.rchord, self.tchord, self.bspan, self.afLoc)
            self.wing.autoWing(15, 15)
            self.aero.set_wing(self.wing)
            self.assertTrue(isinstance(self.aero.surfaces["wing"],Surface))

        def test_build_geometry(self):
            self.aero = Aero()
            self.wing = Wing()
            self.wing.autoSections(self.rchord, self.tchord, self.bspan, self.afLoc)
            self.wing.autoWing(15, 15)
            with self.assertRaises(ValueError):
                self.aero.build_geometry(1,1,1)
            self.aero.set_wing(self.wing)
            self.aero.build_geometry(1,1,1)
            self.assertTrue(isinstance(self.aero.geometry,Geometry))

        def test_add_case(self):
            self.aero = Aero()
            self.wing = Wing()
            with self.assertRaises(ValueError):
                self.aero.add_case()
            self.aero.add_case(name="This is a case",alpha=+2)
            self.aero.add_case(case=Case(name="This is a case",alpha=+2))

        def test_run(self):
            self.aero = Aero()
            self.wing = Wing()
            self.wing.autoSections(self.rchord, self.tchord, self.bspan, self.afLoc)
            self.wing.autoWing(15, 15)
            self.aero.set_wing(self.wing)
            self.aero.build_geometry(1,1,1)
            self.aero.add_case(name="This is a case",alpha=+2)
            self.aero.run()

    unittest.main(verbosity=2)
