import unittest
import modelbase as mb
import numpy as np

class ParameterTests(unittest.TestCase):
    """Tests for parameter initialization and updates"""

    # Initialization
    def test_type_error_parameter_init_list_empty(self):
        m = mb.Model()
        self.assertTrue(isinstance(m.par, mb.parameters.ParameterSet))

    def test_type_error_parameter_init_list(self):
        with self.assertRaises(TypeError):
            m = mb.Model([])

    def test_type_error_parameter_init_int(self):
        with self.assertRaises(TypeError):
            m = mb.Model(1)

    def test_type_error_parameter_init_float(self):
        with self.assertRaises(TypeError):
            m = mb.Model(1.0)

    def test_type_error_parameter_init_str(self):
        with self.assertRaises(TypeError):
            m = mb.Model("a")

    def test_type_error_parameter_init_set(self):
        with self.assertRaises(TypeError):
            m = mb.Model({"a"})

    def test_dictionary_parameter_init(self):
        m = mb.Model()
        p = {"p1": 1}
        m = mb.Model(p)
        self.assertEqual(m.par.p1, p["p1"])

    def test_parameter_set_parameter_init(self):
        m = mb.Model()
        p = mb.parameters.ParameterSet({"p1":1})
        m = mb.Model(p)
        self.assertEqual(m.par.p1, p.p1)

    # Updates
    def test_type_error_parameter_update_list_empty(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.par.update(None)
        with self.assertRaises(TypeError):
            m.par.update()

    def test_type_error_parameter_update_list(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.par.update([])

    def test_type_error_parameter_update_int(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.par.update(1)

    def test_type_error_parameter_update_float(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.par.update(1.0)

    def test_type_error_parameter_update_str(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.par.update("a")

    def test_type_error_parameter_update_set(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.par.update({"a"})

    def test_dictionary_parameter_update(self):
        p = {"p1":1}
        p_update = {"p2":2}
        m = mb.Model(p)
        m.par.update(p_update)
        self.assertEqual(m.par.p1, p["p1"])
        self.assertEqual(m.par.p2, p_update["p2"])

    def test_dictionary_parameter_overwrite(self):
        p = {"p1":1}
        p_update = {"p1":2}
        m = mb.Model(p)
        m.par.update(p_update)
        self.assertEqual(m.par.p1, p_update["p1"])

    def test_parameter_set_parameter_update(self):
        p = {"p1":1}
        p_update = mb.parameters.ParameterSet({"p2":2})
        m = mb.Model(p)
        m.par.update(p_update)
        self.assertEqual(m.par.p1, p["p1"])
        self.assertEqual(m.par.p2, p_update.p2)

    def test_parameter_set_parameter_overwrite(self):
        p = {"p1":1}
        p_update = mb.parameters.ParameterSet({"p1":2})
        m = mb.Model(p)
        m.par.update(p_update)
        self.assertEqual(m.par.p1, p_update.p1)

        # Accesssion
    def test_nonexistent_parameter_accession(self):
        m = mb.Model()
        with self.assertRaises(AttributeError):
            m.par.p1


class CompoundTests(unittest.TestCase):
    """Ensures, that compounds are properly added to the model"""

    def test_set_cpds_type_error_none(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.set_cpds(None)

    def test_set_cpds_type_error_str(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.set_cpds("a")

    def test_set_cpds_type_error_dict(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.set_cpds({})

    def test_set_cpds_type_error_set(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.set_cpds({"a"})

    def test_set_cpds_type_error_int(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.set_cpds(1)

    def test_set_cpds_type_error_float(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.set_cpds(1.0)

    def test_set_cpds_list(self):
        m = mb.Model()
        m.set_cpds([])
        self.assertEqual(m.cpdNames, [])

    def test_set_cpds_list_one_compound(self):
        m = mb.Model()
        m.set_cpds(["X"])
        self.assertEqual(m.cpdNames, ["X"])

    def test_set_cpds_list_two_compounds(self):
        m = mb.Model()
        m.set_cpds(["X", "Y"])
        self.assertEqual(m.cpdNames, ["X", "Y"])

    def test_set_cpds_value_error_double_compounds(self):
        m = mb.Model()
        with self.assertRaises(ValueError):
            m.set_cpds(["X", "X"])

    # Compound ID tests
    def test_compound_ids_empty(self):
        m = mb.Model()
        m.set_cpds([])
        self.assertEqual(m.cpdIdDict, {})
        self.assertEqual(m.cpdIds(), {})

    def test_compound_ids_one_compound(self):
        m = mb.Model()
        m.set_cpds(["X"])
        self.assertEqual(m.cpdIdDict["X"], 0)
        self.assertEqual(m.cpdIds()["X"], 0)

    def test_compound_ids_multiple_compounds(self):
        m = mb.Model()
        m.set_cpds(["X", "Y"])
        self.assertEqual(m.cpdIdDict["X"], 0)
        self.assertEqual(m.cpdIds()["X"], 0)
        self.assertEqual(m.cpdIdDict["Y"], 1)
        self.assertEqual(m.cpdIds()["Y"], 1)

class RatesTests(unittest.TestCase):

    def test_set_rates_rate_name_str(self):
        m = mb.Model()
        m.set_rate("v1", lambda x: x)
        self.assertTrue("v1" in m.rateFn)

    def test_set_rates_rate_name_int(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.set_rate(1, lambda x: x)

    def test_set_rates_rate_name_float(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.set_rate(1.0, lambda x: x)

    def test_set_rates_rate_name_list(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.set_rate([], lambda x: x)

    def test_set_rates_rate_name_dict(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.set_rate({}, lambda x: x)

    def test_set_rates_rate_name_set(self):
        m = mb.Model()
        with self.assertRaises(TypeError):
            m.set_rate({"a"}, lambda x: x)

    # These tests depend on the way we want to construct our
    # rate functions internally. The way we currently do this
    # they do and should not succeed
    # def test_set_rates_func_without_vars(self):
    #     m = mb.Model({"k0":1})
    #     m.set_rate('v0', lambda p: p.k0)
    #     self.assertEqual(m.rateFn["v0"](m.par), m.par.k0)

    # def test_set_rates_func_one_variable(self):
    #     m = mb.Model({"k0":1})
    #     m.set_cpds(["X"])
    #     m.set_rate('v0', lambda p, X: p.k0 * X, "X")
    #     self.assertEqual(m.rateFn["v0"](m.par, 10), m.par.k0 * 10)


class StoichiometryTests(unittest.TestCase):
    """Test basic stoichiometry assembly and stoichiometric matrices"""
    def test_stoichiometries_one_compound_one_reaction(self):
        m = mb.Model({"k":1})
        m.set_cpds(["X"])
        m.set_rate("v0", lambda p: p.k)
        m.set_stoichiometry("v0", {"X":1})
        self.assertEqual(m.stoichiometries, {'v0': {'X': 1}})

    def test_stoichiometries_two_compounds_one_reaction(self):
        m = mb.Model({"k":1})
        m.set_cpds(["X", "Y"])
        m.set_rate("v0", lambda p: p.k)
        m.set_stoichiometry("v0", {"X":-1, "Y":1})
        self.assertEqual(m.stoichiometries, {'v0': {'X': -1, 'Y': 1}})

    def test_stoichiometries_two_compounds_two_reactions(self):
        m = mb.Model({"k":1})
        m.set_cpds(["X", "Y"])
        m.set_rate("v0", lambda p: p.k)
        m.set_rate("v1", lambda p: p.k)
        m.set_stoichiometry("v0", {"X":1})
        m.set_stoichiometry("v1", {"X":-1, "Y":1})
        self.assertEqual(m.stoichiometries, {'v0': {'X': 1}, 'v1': {'X': -1, 'Y': 1}})


class AlgebraicModuleTests(unittest.TestCase):
    """Test basic assembly of algebraic modules"""
    def test_algebraic_module_assembly_one_derived_var(self):
        m = mb.Model()
        m.set_cpds(["X"])
        m.add_algebraicModule(lambda p, x: x, "mod1", ["X"], ["Y"])
        self.assertEqual(m.allCpdNames(), ['X', 'Y'])

    def test_algebraic_module_assembly_two_compounds(self):
        m = mb.Model()
        m.set_cpds(["X", "Y"])
        m.add_algebraicModule(lambda p, x, y: x * y, "mod1", ["X", "Y"], ["Z"])
        self.assertEqual(m.allCpdNames(), ["X", "Y", "Z"])

    def test_algebraic_module_assembly_two_derived_compounds(self):
        m = mb.Model()
        m.set_cpds(["X"])
        m.add_algebraicModule(lambda p, x: (x, x), "mod1", ["X"], ["Y", "Z"])
        self.assertEqual(m.allCpdNames(), ["X", "Y", "Z"])

class RightHandSideTests(unittest.TestCase):
    """Tests functions needed to build the ODE right hand side"""

    def test_full_concentration_vector_type_error_list(self):
        m = mb.Model()
        m.set_cpds(["X"])
        m.add_reaction(
            "v0",
            lambda p, x: x,
            {"X":1},
            "X")
        m.add_algebraicModule(
            lambda p, x: x,
            "mod1",
            ["X"],
            ["Y"])
        with self.assertRaises(TypeError):
            m.fullConcVec([1])

    def test_full_concentration_vector_one_algebraic_module(self):
        m = mb.Model()
        m.set_cpds(["X"])
        m.add_reaction(
            "v0",
            lambda p, x: x,
            {"X":1},
            "X")
        m.add_algebraicModule(
            lambda p, x: x,
            "mod1",
            ["X"],
            ["Y"])
        self.assertEqual(m.fullConcVec(np.array([1])), np.array([1, 1]))

    def test_full_concentration_vector_two_variables(self):
        m = mb.Model()
        m.set_cpds(["X", "Y"])
        m.add_reaction(
            "v0",
            lambda p, x: x,
            {"X":1},
            "X")
        m.add_algebraicModule(
            lambda p, x: x[0] * x[1],
            "mod1",
            ["X", "Y"],
            ["Z"])
        self.assertEqual(m.fullConcVec(np.array([1, 1])), np.array([1, 1, 1]))

    def test_full_concentration_vector_one_two_derived_variables_list(self):
        m = mb.Model()
        m.set_cpds(["X"])
        m.add_reaction(
            "v0",
            lambda p, x: x,
            {"X":1},
            "X")
        m.add_algebraicModule(
            lambda p, x: [x, x],
            "mod1",
            ["X"],
            ["Y", "Z"])
        self.assertEqual(m.fullConcVec(np.array([1])), np.array([1, 1, 1]))

    def test_full_concentration_vector_one_two_derived_variables_tuple(self):
        m = mb.Model()
        m.set_cpds(["X"])
        m.add_reaction(
            "v0",
            lambda p, x: x,
            {"X":1},
            "X")
        m.add_algebraicModule(
            lambda p, x: (x, x),
            "mod1",
            ["X"],
            ["Y", "Z"])
        self.assertEqual(m.fullConcVec(np.array([1])), np.array([1, 1, 1]))

class LabelModelTest(unittest.TestCase):
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

class IntegrationTests(unittest.TestCase):
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
        self.assertTrue(np.isclose(y[-1][0], 22028.621954073271))

        # Two variables
        m = mb.Model({"k0":1, "k1":2})
        m.set_cpds(["X", "Y"])
        m.set_rate('v0', lambda p, x: p.k0 * x, "X")
        m.set_rate('v1', lambda p, x: p.k1 * x, "Y")
        m.set_stoichiometry("v0", {"X":1})
        m.set_stoichiometry("v1", {"Y":1})
        s = mb.Simulate(m)
        y = s.timeCourse(np.linspace(0,10,10),np.ones(2))
        self.assertTrue(np.isclose(y[-1][0], 22026.544796594942))
        self.assertTrue(np.isclose(y[-1][1], 485209580.83282995))

if "__main__" == __name__:
    unittest.main()
