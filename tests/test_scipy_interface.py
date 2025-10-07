from flume.base_classes.system import System
from flume.interfaces.scipy_interface import FlumeScipyInterface
from icecream import ic
import numpy as np
import unittest
from tests.opt_problem_classes import Rosenbrock, RosenbrockConstraint, RosenbrockDVs


class TestUnconstrainedRosenbrock(unittest.TestCase):
    """
    Tests the implementation of the optimization of the unconstrained Rosenbrock function using the FlumeScipyInterface.
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

    def test_optimization_scipy(self):
        """
        Tests the implementation of the SciPy interface for the unconstrained Rosenbrock problem.
        """

        # Construct the interface between Flume and SciPy
        interface = FlumeScipyInterface(flume_sys=self.flume_sys)

        # Set a random starting point
        x0 = np.random.uniform(low=-5.0, high=5.0, size=2)

        # Optimize the problem with SciPy minimize
        x, res = interface.optimize_system(x0=x0)

        # Set the expected optimal values
        xstar = np.array([1.0, 1.0])

        # Check that the values match
        np.testing.assert_allclose(
            actual=x,
            desired=xstar,
            rtol=1e-3,
            err_msg="The computed optimal values do not match the expected solution for the constrained Rosenbrock function.",
            verbose=True,
        )

        return


class TestConstrainedRosenbrock(unittest.TestCase):
    """
    Tests the implementation of the optimization of the constrained Rosenbrock function using the FlumeScipyInterface.
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
            global_con_name={"con.g": {"direction": "leq", "rhs": 1.0}}
        )

        # Store the system as an attribute
        self.flume_sys = sys

        return

    def test_optimization_SLSQP(self):
        """
        Tests that the solution for the constrained Rosenbrock optimization matches the expected solution.
        """

        # Construct the Scipy interface
        interface = FlumeScipyInterface(flume_sys=self.flume_sys)

        # Set a random starting point
        x0 = np.random.uniform(low=-5.0, high=5.0, size=2)

        # Optimize the problem with SciPy minimize
        x, res = interface.optimize_system(x0=x0, method="SLSQP")

        # Set the expected optimal values
        xstar = np.array([0.786, 0.618])

        # Check that the values match
        np.testing.assert_allclose(
            actual=x,
            desired=xstar,
            rtol=1e-3,
            err_msg="The computed optimal values do not match the expected solution for the constrained Rosenbrock function.",
            verbose=True,
        )

        return
