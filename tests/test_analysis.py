import numpy as np
from base_classes.analysis import Analysis
from base_classes.state import State
import unittest

class FunctionAnalysis(Analysis):
    """
    Analysis class that computes the scalar value of an analytic function. This class is used for testing purposes.
    """

    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):
        """
        Analysis class for computing the analytic function for the unittesting.

        f(x, y) = a * x^2 + b * cos(y) + sin(x + y)

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
        self.default_parameters = {"a":1.0, "b":2.0}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Construct the State objects for the variables
        x_var = State(value=2.0, desc="x variable input", source=self)

        y_var = State(value=1.0, desc="y variable input", source=self)

        # Construct the variables dictionary
        self.variables = {"x":x_var, "y":y_var}

        return
    
    def _analyze(self):
        """
        Private analysis method for the object that computes all of the outputs with the provided input variables and parameters to the object
        """

        # Extract the variable and parameter values
        x = self.variables["x"].value
        y = self.variables["y"].value
        a = self.parameters["a"]
        b = self.parameters["b"]

        # Compute the value of f given the variables and parameters
        self.f = a * x**2 + b * np.cos(y) + np.sin(x + y)

        # Update the attribute to reflect that the object has been analyzed
        self.analyzed = True

        # Store the outputs in the output dictionary
        self.outputs = {}

        self.outputs["f"] = State(value=self.f, desc="value of the output", source=self)

        return
    
    def _analyze_adjoint(self):
        """
        Private adjoint analysis method for the object that computes the derivative values for the input variables via the adjoint method.
        """

        # Extract the existing derivative values for the outputs and variables
        fb = self.outputs["f"].deriv

        xb = self.variables["x"].deriv
        yb = self.variables["y"].deriv

        # Extract the variable and parameter values
        x = self.variables["x"].value
        y = self.variables["y"].value
        a = self.parameters["a"]
        b = self.parameters["b"]

        # Compute the contributions to x from f
        xb += fb * (2 * a * x + np.cos(x + y))

        # Compute the contribuitons to y from f
        yb += fb * (-b * np.sin(y) + np.cos(x + y))

        # Set the derivative values for the variables
        self.variables["x"].set_deriv_value(deriv_val=xb)
        self.variables["y"].set_deriv_value(deriv_val=yb)

        # Update the attribute for the adjoint analysis
        self.adjoint_analyzed = True

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

        self.test_obj = FunctionAnalysis(obj_name="test", a=self.a, b=self.b)

        return 
    
    def test_parameter_values(self):
        """
        Checks that the stored parameter values match those randomly computed
        """

        # Extract the parameter names
        param_names = list(self.test_obj.parameters.keys())

        # Store the known parameter values
        param_vals = [self.a, self.b]

        # Loop throuch each parameter and perform the check
        for i, name in enumerate(param_names):
            with self.subTest(f"parameter = {name}"):
                self.assertEqual(param_vals[i], self.test_obj.parameters[name], f"Values for parameter {name} do not match")

        return
    
    def test_set_var_values(self):
        """
        Checks the functionality of the set variable values method
        """

        # Set the variable values for the analysis object
        self.test_obj.set_var_values(variables={"x":self.x, "y":self.y})

        # Store the known variable values
        var_names = ["x", "y"]
        known_vals = {"x":self.x, "y":self.y}

        # Extract the variable values
        var_vals = self.test_obj.get_var_values()

        # Loop through each varaible and perform the check
        for i, name in enumerate(var_names):
            with self.subTest(f"variable = {name}"):
                self.assertEqual(known_vals[name], var_vals[name], f"Values for variable {name} do not match")

        return
    
    def test_output_values(self):
        """
        Checks the computation of the output value for the analytic function
        """

        # Set the variable values for the analysis object
        self.test_obj.set_var_values(variables={"x":self.x, "y":self.y})

        # Analyze the system
        self.test_obj.analyze()

        # Extract the output value
        outputs = self.test_obj.get_output_values(outputs=["f"])
        f_obj = outputs["f"]

        # Compute f with the known inputs
        f_exp = self.a * self.x**2 + self.b * np.cos(self.y) + np.sin(self.x + self.y)

        # Check the values
        self.assertEqual(f_obj, f_exp, "Output values for f do not match")

        return
    
    def test_isolated_adjoint(self):
        """
        Checks the computation of the isolated adjoint method for the analysis object
        """

        # Set the variable values for the analysis object
        self.test_obj.set_var_values(variables={"x":self.x, "y":self.y})

        # Declare both variables as design variables
        self.test_obj.declare_design_vars(variables=["x", "y"])

        # Perform the combined adjoint test for the object
        rel_error = self.test_obj.test_isolated_adjoint(print_res=False, debug_print=False, method="cs", defined_vars={"test.x":2.0})

        # Check to make sure the relative error is less than 1e-6
        self.assertLessEqual(rel_error, 1e-5, "Relative error for the isolated adjoint test is greater than 1e-5.")

        return    
        
    def test_combined_adjoint(self):
        """
        Checks the computation of the combined adjoint method for the analysis object
        """

        # Set the variable values for the analysis object
        self.test_obj.set_var_values(variables={"x":self.x, "y":self.y})

        # Declare both variables as design variables
        self.test_obj.declare_design_vars(variables=["x", "y"])

        # Perform the combined adjoint test for the object
        rel_error = self.test_obj.test_combined_adjoint(print_res=False, debug_print=False, method="cs", defined_vars={"test.x":2.0})

        # Check to make sure the relative error is less than 1e-6
        self.assertLessEqual(rel_error, 1e-5, "Relative error for the combined adjoint test is greater than 1e-5.")

        return    

if __name__ == "__main__":
    # Perform the unittests
    unittest.main()
