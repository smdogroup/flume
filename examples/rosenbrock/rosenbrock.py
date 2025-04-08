from flume.base_classes.analysis import Analysis
from flume.base_classes.state import State
from flume.base_classes.system import System
from flume.interfaces.paropt_interface import FlumeParOptInterface
from icecream import ic
from paropt import ParOpt
import numpy as np


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


if __name__ == "__main__":

    # Construct the analysis object for the Rosenbrock function
    a = 1.0
    b = 100.0

    rosenbrock = Rosenbrock(obj_name="rosenbrock", sub_analyses=[], a=a, b=b)

    # Construct the system
    sys = System(
        sys_name="rosen_sys",
        top_level_analysis_list=[rosenbrock],
        log_name="flume.log",
        log_prefix="examples/rosenbrock/test_paropt",
    )

    # Declare the design variables for the system
    sys.declare_design_vars(
        global_var_name={
            "rosenbrock.x": {"lb": -2.0, "ub": 2.0},
            "rosenbrock.y": {"lb": -2.0, "ub": 2.0},
        }
    )

    sys.declare_objective(global_obj_name="rosenbrock.f")

    # Create the FlumeParOptInterface
    interface = FlumeParOptInterface(flume_sys=sys)

    # Construct the paropt problem for the Flume system
    paroptprob = interface.construct_paropt_problem()
    options = interface.get_paropt_default_options(
        output_prefix="examples/rosenbrock/test_paropt"
    )

    # Perform the optimization with ParOpt
    opt = ParOpt.Optimizer(paroptprob, options)
    opt.optimize()

    # Extract the optimized point
    x, z, zw, zl, zu = opt.getOptimizedPoint()
    x = np.array(x)

    # Check that the values of the optimized point match expected values
    xstar = 1.0
    ystar = 1.0

    # Compute the rel error
    xerr = (xstar - x[0]) / xstar
    yerr = (ystar - x[1]) / ystar

    print("\n%15s %15s %15s %15s" % ("Variable", "Expected", "Computed", "Error"))
    print("\n%15s %15s %15s %15s" % ("x", xstar, x[0], xerr))
    print("\n%15s %15s %15s %15s" % ("y", ystar, x[1], yerr))
