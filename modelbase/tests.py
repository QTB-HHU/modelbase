import unittest
import modelbase as mb
import numpy as np

class AssemblyTests(unittest.TestCase):
    """Testing user input for model assembly"""

    def test_parameter_init(self):
        """
        Ensures that the model can be initialized properly
        """
        # No parameters
        m = mb.Model()
        # how can we test this?

        # Dictionary
        p = {"p1": 1}
        m = mb.Model(p)
        assert m.par.p1 == 1

        # ParameterSet
        p = mb.parameters.ParameterSet({"p1":1})
        m = mb.Model(p)
        assert m.par.p1 == 1


    def test_set_cpds(self):
        """
        Ensures, that compounds are properly added
        to the model
        """
        # No compounds
        m = mb.Model()
        m.set_cpds([])
        assert m.cpdNames == []

        # One compound
        m = mb.Model()
        m.set_cpds(["X"])
        assert m.cpdNames == ["X"]

        # Multiple compounds
        m = mb.Model()
        m.set_cpds(["X", "Y"])
        assert m.cpdNames == ["X", "Y"]

        # Double compounds
        m = mb.Model()
        m.set_cpds(["X", "X"])
        assert m.cpdNames == ["X"]

    def test_compound_ids(self):
        # No compounds
        m = mb.Model()
        m.set_cpds([])
        assert m.cpdIdDict == {}
        assert m.cpdIds() == {}

        # One compound
        m = mb.Model()
        m.set_cpds(["X"])
        assert m.cpdIdDict["X"] == 0
        assert m.cpdIds()["X"] == 0

        # Multiple compounds
        m = mb.Model()
        m.set_cpds(["X", "Y"])
        assert m.cpdIdDict["X"] == 0
        assert m.cpdIdDict["Y"] == 1
        assert m.cpdIds()["X"] == 0
        assert m.cpdIds()["Y"] == 1

        # Double compounds
        m = mb.Model()
        m.set_cpds(["X", "X"])
        assert m.cpdIdDict["X"] == 0
        assert m.cpdIds()["X"] == 0

    def test_set_rates(self):
        # Reactions without compounds
        m = mb.Model({"k0":1})
        m.set_rate('v0',lambda p: p.k0)
        assert m.rateFn["v0"](m.par) == 1

        ## Reactions with multiple compounds, do we want to have this?
        # m = mb.Model({"k0":10})
        # m.set_rate('v0',lambda p, x: p.k0 * x)
        # assert m.rateFn["v0"](m.par, 10) == 100
        m = mb.Model({"k0":10})
        m.set_cpds(["X"])
        m.set_rate('v0', lambda p, x: p.k0*x, "X")
        m.set_stoichiometry("v0", {"X":1})
        assert m.rates(np.array([10]))["v0"] == 100

    def test_stoichiometries(self):
        """
        Test basic stoichiometry assembly and stoichiometric matrices
        """
        # One compound
        m = mb.Model({"k":1})
        m.set_cpds(["X"])
        m.set_rate("v0", lambda p: p.k)
        m.set_stoichiometry("v0", {"X":1})
        assert m.stoichiometries == {'v0': {'X': 1}}
        # Is there a way to test the matrices?
        # assert m.stoichiometryMatrix() == np.matrix([[ 1.]])

        # Two compounds, two reactions
        m = mb.Model({"k":1})
        m.set_cpds(["X", "Y"])
        m.set_rate("v0", lambda p: p.k)
        m.set_rate("v1", lambda p: p.k)
        m.set_stoichiometry("v0", {"X":1})
        m.set_stoichiometry("v1", {"X":-1, "Y":1})
        assert m.stoichiometries == {'v0': {'X': 1}, 'v1': {'X': -1, 'Y': 1}}
        # Is there a way to test the matrices?
        # assert m.stoichiometryMatrix() == np.matrix([[ 1., -1.],[ 0.,  1.]])

    def test_algebraic_module_assembly(self):
        m = mb.Model()
        m.set_cpds(["X"])
        m.add_algebraicModule(lambda x: x, "mod1", ["X"], ["Y"])
        assert m.allCpdNames() == ['X', 'Y']

    def test_integration(self):
        """
        Test basic integration output
        """
        # One variable
        m = mb.Model({"k0":1})
        m.set_cpds(["X"])
        m.set_rate('v0', lambda p, x: p.k0 * x, "X")
        m.set_stoichiometry("v0", {"X":1})
        s = mb.Simulate(m)
        y = s.timeCourse(np.linspace(0,10,10),np.ones(1))
        assert y[-1][0] == 22028.621954073271

        # Two variables
        m = mb.Model({"k0":1, "k1":2})
        m.set_cpds(["X", "Y"])
        m.set_rate('v0', lambda p, x: p.k0 * x, "X")
        m.set_rate('v1', lambda p, x: p.k1 * x, "Y")
        m.set_stoichiometry("v0", {"X":1})
        m.set_stoichiometry("v1", {"Y":1})
        s = mb.Simulate(m)
        y = s.timeCourse(np.linspace(0,10,10),np.ones(2))
        assert y[-1][0] == 22026.544796594942
        assert y[-1][1] == 485209580.83282995

    def test_label_compound_assembly(self):
        m = mb.LabelModel()
        m.add_base_cpd("A", 2)
        assert m.allCpdNames() == ['A00', 'A01', 'A10', 'A11', 'A_total']

    def test_carbonmap_reaction_assembly(self):
        m = mb.LabelModel()
        m.add_base_cpd("A", 2)
        m.add_base_cpd("B", 2)
        m.add_carbonmap_reaction("v1", lambda x: x, [0, 1], ["A"], ["B"], "A")
        m.add_carbonmap_reaction("v2", lambda x: x, [0, 1], ["B"], ["A"], "B")
        list(m.rateFn.keys()) == ['v100', 'v101', 'v110', 'v111',
                                  'v200', 'v201', 'v210', 'v211']


if "__main__" == __name__:
    unittest.main()
