import numpy as np
from flume.base_classes.analysis import Analysis
from flume.base_classes.state import State
import unittest
from jax import numpy as jnp

class Aobj1(Analysis):
    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):
        # Set the default parameters for the object
        self.default_parameters = {"a": 1.0, "b": 2.0}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Construct the State objects for the variables
        x_var = State(value=2.0, desc="x variable input", source=self)

        # Construct the variables dictionary
        self.variables = {"x": x_var}
        return

    def callable_analyze(self,x:list,a:list):
        f = x[0]**2
        self.ca_descs = {"f":"function 1 value"}
        return f

class Aobj2(Analysis):
    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):
        # Set the default parameters for the object
        self.default_parameters = {"a": 1.0, "b": 2.0}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Construct the State objects for the variables
        f_var = State(value=2.0, desc="x variable input", source=self)
        """ HERE """

        # Construct the variables dictionary
        self.variables = {"f": f_var}
        return

    def callable_analyze(self,x:list,a:list):
        g = jnp.log(x[0])
        self.ca_descs = {"g":"function 2 value"}
        return g

class TestAnalysisObject(unittest.TestCase):
    """
    Tests the implementation of an Analysis object for a simple analytic function.
    """

    def setUp(self):
        """
        Sets up the sample Analysis object that is used for the subsequent tests.
        """

        # Set random values for the parameter and state values
        self.a = np.random.uniform()
        self.b = np.random.uniform()

        self.x = np.random.uniform(low=0.0, high=np.pi)

        self.obj_1 = Aobj1("obj1",sub_analyses = [])
        self.obj_2 = Aobj2("obj2",sub_analyses = [self.obj_1])

        return

    def test_subanalyses(self):
        self.obj_1.declare_design_vars(variables=["x"])
        adj_err = self.obj_2.test_combined_adjoint(debug_print=True)

        # Check to make sure the relative error is less than 1e-5
        self.assertLessEqual(
            adj_err,
            1e-5,
            "Relative error for the adjoint test is greater than 1e-5.",
        )

if __name__ == "__main__":
    # Perform the unittests
    unittest.main()
