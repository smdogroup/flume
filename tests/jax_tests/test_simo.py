import numpy as np
from flume.base_classes.analysis import Analysis
from flume.base_classes.state import State
import unittest
from jax import numpy as jnp

class FunctionAnalysis_analysis(Analysis):
    """
    Analysis class that computes the scalar value of an analytic function. This class is used for testing purposes.
    """
    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):
        """
        Analysis class for computing the analytic function for the unittesting.


        Parameters
        ----------
        obj_name : str
            Name for the analysis object
        sub_analyses : list
            A list of any sub-analyses for the DeltaV analysis (nominally, this is empty)

        Keyword Arguments
        -----------------
        a : float
            First parameter for the analytic function
        b : float
            Second parameter for the analytic function
        """

        # Set the default parameters for the object
        self.default_parameters = {"a": 1.0, "b": 2.0}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Construct the State objects for the variables
        x_var = State(value=2.0, desc="x variable input", source=self)
        y_var = State(value=1.0, desc="y variable input", source=self)

        # Construct the variables dictionary
        self.variables = {"x": x_var, "y": y_var}

        return

    def _analyze(self):
        """
        Private analysis method for the object that computes all of the outputs with the provided input variables and parameters to the object
        """

        # Extract the variable and parameter values
        x = self.variables["x"].value
        a = self.parameters["a"]
        b = self.parameters["b"]

        # Compute the value of f given the variables and parameters
        self.f = np.sin(x)*a*b**2
        self.g = a*x**2-b

        # Update the attribute to reflect that the object has been analyzed
        self.analyzed = True

        # Store the outputs in the output dictionary
        self.outputs = {}

        self.outputs["f"] = State(value=self.f, desc="value of output 1", source=self)
        self.outputs["g"] = State(value=self.f, desc="value of output 2", source=self)

        return

    def _analyze_adjoint(self):
        """
        Private adjoint analysis method for the object that computes the derivative values for the input variables via the adjoint method.
        """

        # Extract the existing derivative values for the outputs and variables
        fb = self.outputs["f"].deriv
        gb = self.outputs["g"].deriv

        xb = self.variables["x"].deriv
        yb = self.variables["y"].deriv

        # Extract the variable and parameter values
        x = self.variables["x"].value
        y = self.variables["y"].value
        a = self.parameters["a"]
        b = self.parameters["b"]

        # Compute the contributions to x from f
        xb += fb * np.cos(x)*a*b**2
        # Compute the contribuitons to y from f
        yb += 0

        # Compute the contributions to x from f
        xb += gb * 2*a*x
        # Compute the contribuitons to y from f
        yb += gb * -1

        # Set the derivative values for the variables
        self.variables["x"].set_deriv_value(deriv_val=xb)
        self.variables["y"].set_deriv_value(deriv_val=yb)

        # Update the attribute for the adjoint analysis
        self.adjoint_analyzed = True

class FunctionAnalysis_wrapper(Analysis):
    """
    Analysis class that computes the scalar value of an analytic function. This class is used for testing purposes.
    """

    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):
        """
        Analysis class for computing the analytic function for the unittesting.

        Parameters
        ----------
        obj_name : str
            Name for the analysis object
        sub_analyses : list
            A list of any sub-analyses for the DeltaV analysis (nominally, this is empty)

        Keyword Arguments
        -----------------
        a : float
            First parameter for the analytic function
        b : float
            Second parameter for the analytic function
        """

        # Set the default parameters for the object
        self.default_parameters = {"a": 1.0, "b": 2.0}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Construct the State objects for the variables
        x_var = State(value=2.0, desc="x variable input", source=self)
        y_var = State(value=1.0, desc="y variable input", source=self)

        # Construct the variables dictionary
        self.variables = {"x": x_var, "y": y_var}

        return

    def callable_analyze(self,x: list,a: list):
        """
        Private analysis method for the object that computes all of the outputs with the provided input variables and parameters to the object
        """
        # Compute the value of f given the variables and parameters
        f = jnp.sin(x[0])*a[0]*a[1]**2
        g = 2*a[0]*x[0]-a[1]
        self.ca_descs = {"f":"function 1 value","g":"function 2 value"}
        return f,g

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
        self.y = np.random.uniform(low=0.0, high=np.pi)

        self.test_obj_1 = FunctionAnalysis_analysis(obj_name="test", a=self.a, b=self.b)
        self.test_obj_2 = FunctionAnalysis_wrapper(obj_name="test", a=self.a, b=self.b)

        return

    def test_parameter_values(self):
        """
        Checks that the stored parameter values match those randomly computed
        """

        # Extract the parameter names
        param_names_1 = self.test_obj_1.parameters
        param_names_2 = self.test_obj_2.parameters

        # Loop throuch each parameter and perform the check
        for i, (name1,name2) in enumerate(zip(param_names_1,param_names_2)):
            with self.subTest(f"parameter = {i}, name check:"):
                self.assertEqual(name1,name2,
                                 f"Name {name1} and name {name2} do not match")
            with self.subTest(f"parameter = {i}, value check:"):
                self.assertEqual(param_names_1[name1],param_names_2[name2],
                                 f"Analyze value and callable analyze value do not match")

        return

    def test_set_var_values(self):
        """
        Checks the functionality of the set variable values method
        """

        # Set the variable values for the analysis object
        self.test_obj_1.set_var_values(variables={"x": self.x, "y": self.y})
        self.test_obj_2.set_var_values(variables={"x": self.x, "y": self.y})

        # Extract the variable values
        var_vals_1 = self.test_obj_1.get_var_values()
        var_vals_2 = self.test_obj_1.get_var_values()

        # Loop through each varaible and perform the check
        for i, name in enumerate(var_vals_1):
            with self.subTest(f"variable = {name}"):
                self.assertEqual(
                    var_vals_1[name],
                    var_vals_2[name],
                    f"Values for variable {name} do not match",
                )
        return

    def test_output_values(self):
        """
        Checks the computation of the output value for the analytic function
        """

        # Set the variable values for the analysis object
        self.test_obj_1.set_var_values(variables={"x": self.x, "y": self.y})
        self.test_obj_2.set_var_values(variables={"x": self.x, "y": self.y})

        # Analyze the system
        self.test_obj_1.analyze()
        self.test_obj_2.analyze()

        # Extract the output value
        outputs1 = self.test_obj_1.get_output_values(outputs=["f"])
        f_obj_1 = outputs1["f"]
        outputs2 = self.test_obj_2.get_output_values(outputs=["f"])
        f_obj_2 = outputs2["f"]

        # Check the values
        self.assertAlmostEqual(f_obj_1, f_obj_2, places=5, msg="Output values for f do not match")

        return

    def test_adjoint_values(self):
        """
        Checks the computation of the output value for the analytic function
        """
        self.test_obj_2.declare_design_vars(variables=["x","y"])
        rel_error = self.test_obj_2.test_combined_adjoint()

        # Check to make sure the relative error is less than 1e-5
        self.assertLessEqual(
            rel_error,
            1.e-5,
            "Relative error for the isolated adjoint test is greater than 1e-5.",
        )

        return

if __name__ == "__main__":
    # Perform the unittests
    unittest.main()
