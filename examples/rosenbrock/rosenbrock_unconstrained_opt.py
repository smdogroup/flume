# Description: this defines the run script that combines the individual Analysis classes to perform optimization for the unconstrained Rosenbrock problem.


from flume.base_classes.system import System
from icecream import ic
import numpy as np
import argparse
from .rosenbrock_problem_classes import Rosenbrock


class RosenbrockCallback:
    def __init__(self, extra_param, obj=None):
        self.extra_param = extra_param
        self.obj = obj

    def __call__(self, x, it_num):
        print(f"Iteration {it_num}: x, y = {x[0], x[1]}")


if __name__ == "__main__":

    # Construct the argument parser
    parser = argparse.ArgumentParser(
        description="Program performs the optimization of the unconstrained Rosenbrock function, using the FlumeScipyInterface by default."
    )

    # Add an argument that specifies whether ParOpt should be used to perform the optimization, otherwise SciPy is used.
    parser.add_argument(
        "--use-paropt",
        "-paropt",
        action="store_true",
        help="Argument that, when provided, triggers the script to utilize ParOpt to perform the optimization. Otherwise, SciPy is used.",
    )

    args = parser.parse_args()

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
        log_prefix="examples/rosenbrock/unconstrained",
    )

    # Declare the design variables for the system
    sys.declare_design_vars(
        global_var_name={
            "rosenbrock.x": {"lb": -2.0, "ub": 2.0},
            "rosenbrock.y": {"lb": -2.0, "ub": 2.0},
        }
    )

    # Declare the objective function
    sys.declare_objective(global_obj_name="rosenbrock.f")

    # Perform the optimization, using either SciPy or ParOpt depending on the input arguments
    if args.use_paropt:
        from flume.interfaces.paropt_interface import FlumeParOptInterface
        from paropt import ParOpt

        # Construct the callback function to use
        callback = RosenbrockCallback(extra_param="test", obj=sys)

        # Create the FlumeParOptInterface
        interface = FlumeParOptInterface(flume_sys=sys, callback=callback)

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
    else:
        from flume.interfaces.scipy_interface import FlumeScipyInterface

        # Construct the FlumeScipyInterface
        interface = FlumeScipyInterface(flume_sys=sys)

        # Set the initial point
        x0 = np.array([0.0, 0.0])

        # Perform the optimization
        x, res = interface.optimize_system(x0=x0, method="SLSQP", maxit=100)

    x = np.array(x)

    # Check that the values of the optimized point match expected values
    xstar = 1.0
    ystar = 1.0

    # Compute the rel error
    xerr = (xstar - x[0]) / xstar
    yerr = (ystar - x[1]) / ystar

    print("\n%15s %15s %15s %15s" % ("Variable", "Expected", "Computed", "Error"))
    print("%15s %15s %15s %15s" % ("x", xstar, x[0], xerr))
    print("%15s %15s %15s %15s" % ("y", ystar, x[1], yerr))
