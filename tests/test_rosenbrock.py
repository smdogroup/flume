from flume.base_classes.analysis import Analysis
from flume.base_classes.state import State
from flume.base_classes.system import System
from flume.interfaces.paropt_interface import FlumeParOptInterface
from icecream import ic
from paropt import ParOpt
import numpy as np
import unittest


class Rosenbrock(Analysis):

    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):

        # Set the default parameters
        self.default_parameters = {"a": 1.0, "b": 100.0}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the default State for the variables
        xvar = State(value=0.0, desc="x state value", source=self)
        yvar = State(value=0.0, desc="y state value", source=self)

        # Construct variables dictionary
        self.variables = {"x": xvar, "y": yvar}

        return

    def _analyze(self):
        # Extract the variable values
        x = self.variables["x"].value
        y = self.variables["y"].value

        # Extract the parameter values
        a = self.parameters["a"]
        b = self.parameters["b"]

        # Compute the value of the Rosenbrock function
        f = (a - x) ** 2 + b * (y - x**2) ** 2

        # Update the analyzed attribute
        self.analyzed = True

        # Store the outputs
        self.outputs = {}

        self.outputs["f"] = State(
            value=f, desc="Rosenbrock function value", source=self
        )

        return

    def _analyze_adjoint(self):
        # Extract the derivatives of the outputs
        fb = self.outputs["f"].deriv

        # Extract the variable values
        x = self.variables["x"].value
        y = self.variables["y"].value

        # Extract the parameter values
        a = self.parameters["a"]
        b = self.parameters["b"]

        # Compute xb
        xb = 2 * (a - x) * -1 + 2.0 * b * (y - x**2) * -2 * x

        # Compute yb
        yb = 2 * b * (y - x**2)

        # Update the analyzed adjoint attribute
        self.adjoint_analyzed = True

        # Assign the derivative values
        self.variables["x"].set_deriv_value(deriv_val=xb)
        self.variables["y"].set_deriv_value(deriv_val=yb)

        return


class RosenbrockDVs(Analysis):
    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):
        # Set the default parameters
        self.default_parameters = {}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the default states for the variables
        x_dv_var = State(
            value=0.0, desc="Value to use for the x design variable", source=self
        )

        y_dv_var = State(
            value=0.0, desc="Value to use for the y design variable", source=self
        )

        self.variables = {"x_dv": x_dv_var, "y_dv": y_dv_var}

        return

    def _analyze(self):

        # Extract the variables
        x_dv = self.variables["x_dv"].value
        y_dv = self.variables["y_dv"].value

        # Update the analyzed attribute
        self.analyzed = True

        # Store the outputs in the outputs dictionary
        self.outputs = {}

        self.outputs["x"] = State(value=x_dv, desc="Value for x", source=self)

        self.outputs["y"] = State(value=y_dv, desc="Value for y", source=self)

        return

    def _analyze_adjoint(self):

        # Extract the output derivatives
        xb = self.outputs["x"].deriv
        yb = self.outputs["y"].deriv

        # Extract the variable derivatives
        x_dvb = self.variables["x_dv"].deriv
        y_dvb = self.variables["y_dv"].deriv

        # Update the derivative values
        x_dvb += xb
        y_dvb += yb

        # Update the analyzed adjoint attribute
        self.adjoint_analyzed = True

        # Set the derivative values
        self.variables["x_dv"].set_deriv_value(deriv_val=x_dvb)

        self.variables["y_dv"].set_deriv_value(deriv_val=y_dvb)

        return


class RosenbrockConstraint(Analysis):
    def __init__(self, obj_name: str, sub_analyses=[RosenbrockDVs], **kwargs):
        # Set the default parameters
        self.default_parameters = {}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the default State for the variables
        xvar = State(value=0.0, desc="x state value", source=self)
        yvar = State(value=0.0, desc="y state value", source=self)

        # Construct variables dictionary
        self.variables = {"x": xvar, "y": yvar}

        return

    def _analyze(self):

        # Extract the variable values
        x = self.variables["x"].value
        y = self.variables["y"].value

        # Compute the value of the constraint (negative sign is used to test constraints that have a negative rhs)
        g = (x**2 + y**2) * -1.0

        # Update the analyzed attribute
        self.analyzed = True

        # Store the outputs in the outputs dictionary
        self.outputs = {}

        self.outputs["g"] = State(
            value=g,
            desc="Constraint value that constraints the design space to a circle",
            source=self,
        )

        return

    def _analyze_adjoint(self):

        # Extract the variable values
        x = self.variables["x"].value
        y = self.variables["y"].value

        # Extract the variable derivatives
        xb = self.variables["x"].deriv
        yb = self.variables["y"].deriv

        # Extract the output derivatives
        gb = self.outputs["g"].deriv

        # Add the contributions to xb and yb
        xb += gb * 2 * x * -1.0
        yb += gb * 2 * y * -1.0

        # Update the analyzed adjoint attribute
        self.adjoint_analyzed = True

        # Assign the derivative values
        self.variables["x"].set_deriv_value(deriv_val=xb)
        self.variables["y"].set_deriv_value(deriv_val=yb)

        return


class TestUnconstrainedRosenbrock(unittest.TestCase):
    """
    Tests the implementation of the optimization of the unconstrained Rosenbrock function using the FlumeParOptInterface.
    """

    def setUp(self):

        # Construct the analysis object for the Rosenbrock function
        a = 1.0
        b = 100.0

        rosenbrock = Rosenbrock(obj_name="rosenbrock", sub_analyses=[], a=a, b=b)

        rosenbrock.set_var_values(variables={"x": 0.0, "y": 0.0})

        # Construct the system
        sys = System(
            sys_name="rosen_sys",
            top_level_analysis_list=[rosenbrock],
            log_name="flume.log",
            log_prefix="tests/rosenbrock_unconstrained",
        )

        # Declare the design variables for the system
        sys.declare_design_vars(
            global_var_name={
                "rosenbrock.x": {"lb": -2.0, "ub": 2.0},
                "rosenbrock.y": {"lb": -2.0, "ub": 2.0},
            }
        )

        sys.declare_objective(global_obj_name="rosenbrock.f")

        # Save the system as an attribute
        self.flume_sys = sys

        return

    def test_optimization(self):
        """
        Tests that the solution for the unconstrained Rosenbrock optimization matches the expected solution.
        """

        interface = FlumeParOptInterface(flume_sys=self.flume_sys)

        # Construct the paropt problem for the Flume system
        paroptprob = interface.construct_paropt_problem()
        options = interface.get_paropt_default_options(
            output_prefix="tests/rosenbrock_unconstrained"
        )

        # Perform the optimization with ParOpt
        opt = ParOpt.Optimizer(paroptprob, options)
        opt.optimize()

        # Extract the optimized point
        x, z, zw, zl, zu = opt.getOptimizedPoint()
        x = np.array(x)

        # Set the expected optimal values
        xstar = np.array([1.0, 1.0])

        # Check that the values match
        np.testing.assert_allclose(
            actual=x,
            desired=xstar,
            rtol=5e-4,
            err_msg="The computed optimal values do not match the expected solution for the unconstrained Rosenbrock function.",
            verbose=True,
        )

        return


class TestConstrainedRosenbrock(unittest.TestCase):
    """
    Tests the implementation of the optimization of the constrained Rosenbrock function using the FlumeParOptInterface.
    """

    def setUp(self):

        # Construct the design variables object
        rosenbrock_dvs = RosenbrockDVs(obj_name="dvs", sub_analyses=[])

        rosenbrock_dvs.set_var_values(variables={"x_dv": 0.0, "y_dv": 0.0})

        # Construct the analysis object for the Rosenbrock function
        a = 1.0
        b = 100.0

        rosenbrock = Rosenbrock(
            obj_name="rosenbrock", sub_analyses=[rosenbrock_dvs], a=a, b=b
        )

        # Construct the analysis object for the constraint on the design variables
        rosenbrock_con = RosenbrockConstraint(
            obj_name="con", sub_analyses=[rosenbrock_dvs]
        )

        # Construct the system
        sys = System(
            sys_name="rosen_sys_con",
            top_level_analysis_list=[rosenbrock, rosenbrock_con],
            log_name="flume.log",
            log_prefix="tests/rosenbrock_constrained",
        )

        # Declare the design variables for the system
        sys.declare_design_vars(
            global_var_name={
                "dvs.x_dv": {"lb": -1.5, "ub": 1.5},
                "dvs.y_dv": {"lb": -1.5, "ub": 1.5},
            }
        )

        # Declare the objective
        sys.declare_objective(global_obj_name="rosenbrock.f")

        # Declare the constraint
        sys.declare_constraints(
            global_con_name={"con.g": {"direction": "geq", "rhs": -2.0}}
        )

        # Store the system as an attribute
        self.flume_sys = sys

        return

    def test_optimization(self):
        """
        Tests that the solution for the constrained Rosenbrock optimization matches the expected solution.
        """

        interface = FlumeParOptInterface(flume_sys=self.flume_sys)

        # Construct the paropt problem for the Flume system
        paroptprob = interface.construct_paropt_problem()
        options = interface.get_paropt_default_options(
            output_prefix="tests/rosenbrock_constrained"
        )

        # Perform the optimization with ParOpt
        opt = ParOpt.Optimizer(paroptprob, options)
        opt.optimize()

        # Extract the optimized point
        x, z, zw, zl, zu = opt.getOptimizedPoint()
        x = np.array(x)

        # Set the expected optimal values
        xstar = np.array([1.0, 1.0])

        # Check that the values match
        np.testing.assert_allclose(
            actual=x,
            desired=xstar,
            rtol=5e-4,
            err_msg="The computed optimal values do not match the expected solution for the constrained Rosenbrock function.",
            verbose=True,
        )

        return


if __name__ == "__main__":
    # Run the unittests
    unittest.main()
