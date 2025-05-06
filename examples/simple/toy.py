from flume.base_classes.analysis import Analysis
from flume.base_classes.state import State
from flume.base_classes.system import System
from flume.interfaces.paropt_interface import FlumeParOptInterface, ParOptProb
from icecream import ic
from paropt import ParOpt
import numpy as np


class Block1(Analysis):

    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):

        # Set default parameters
        self.default_parameters = {}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the default State for the variables
        xvar = State(value=1.0, desc="x state value", source=self)

        # Construct variables dictionary
        self.variables = {"x": xvar}

    def _analyze(self):
        # Extract the variable value
        x = self.variables["x"].value

        # Compute the outputs
        y = 2 * x
        q = x / 2

        # Update analyzed attribute
        self.analyzed = True

        # Store the outputs
        self.outputs = {}

        self.outputs["y"] = State(value=y, desc="y state value", source=self)
        self.outputs["q"] = State(value=q, desc="q state value", source=self)

    def _analyze_adjoint(self):
        # Extract the derivatives of the outputs
        yb = self.outputs["y"].deriv
        qb = self.outputs["q"].deriv

        # Compute the contributions to xb
        xb = yb * 2.0
        xb += qb * 0.5

        # Update the analyzed adjoint attribute
        self.adjoint_analyzed = True

        # Update the derivative values for the varialbes
        self.variables["x"].set_deriv_value(deriv_val=xb)

        return


class Block2(Analysis):

    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):

        # Set default parameters
        self.default_parameters = {}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the default State for the variables
        yvar = State(value=1.0, desc="y state value", source=self)

        # Construct variables dictionary
        self.variables = {"y": yvar}

        return

    def _analyze(self):
        # Extract the variable value
        y = self.variables["y"].value

        # Compute the outputs
        z = y**2

        # Update analyzed attribute
        self.analyzed = True

        # Store the outputs
        self.outputs = {}

        self.outputs["z"] = State(value=z, desc="z state value", source=self)

        return

    def _analyze_adjoint(self):
        # Extract the derivatives of the outputs
        zb = self.outputs["z"].deriv

        # Extract the variable values
        y = self.variables["y"].value

        # Compute yb
        yb = zb * 2.0 * y

        # Update the analyzed adjoint attribute
        self.adjoint_analyzed = True

        # Update the derivative values for the varialbes
        self.variables["y"].set_deriv_value(deriv_val=yb)

        return


class Block3(Analysis):

    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):

        # Set default parameters
        self.default_parameters = {}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the default State for the variables
        zvar = State(value=1.0, desc="z state value", source=self)

        # Construct variables dictionary
        self.variables = {"z": zvar}

        return

    def _analyze(self):
        # Extract the variable value
        z = self.variables["z"].value

        # Compute the outputs
        alpha = z * 2.0

        # Update analyzed attribute
        self.analyzed = True

        # Store the outputs
        self.outputs = {}

        self.outputs["alpha"] = State(
            value=alpha, desc="alpha state value", source=self
        )

        return

    def _analyze_adjoint(self):
        # Extract the derivatives of the outputs
        alphab = self.outputs["alpha"].deriv

        # Compute zb
        zb = alphab * 2.0

        # Update the analyzed adjoint attribute
        self.adjoint_analyzed = True

        # Update the derivative values for the varialbes
        self.variables["z"].set_deriv_value(deriv_val=zb)

        return


class Block4(Analysis):

    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):

        # Set default parameters
        self.default_parameters = {}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the default State for the variables
        qvar = State(value=1.0, desc="q state value", source=self)

        # Construct variables dictionary
        self.variables = {"q": qvar}

        return

    def _analyze(self):
        # Extract the variable value
        q = self.variables["q"].value

        # Compute the outputs
        beta = 1 / q

        # Update analyzed attribute
        self.analyzed = True

        # Store the outputs
        self.outputs = {}

        self.outputs["beta"] = State(value=beta, desc="beta state value", source=self)

        return

    def _analyze_adjoint(self):
        # Extract the derivatives of the outputs
        betab = self.outputs["beta"].deriv

        # Extract the variable values
        q = self.variables["q"].value

        # Compute qb
        qb = betab * -1.0 / q**2

        # Update the analyzed adjoint attribute
        self.adjoint_analyzed = True

        # Update the derivative values for the varialbes
        self.variables["q"].set_deriv_value(deriv_val=qb)

        return


if __name__ == "__main__":

    # Construct the different analysis objects
    B1 = Block1(obj_name="block1")
    B2 = Block2(obj_name="block2", sub_analyses=[B1])
    B3 = Block3(obj_name="block3", sub_analyses=[B2])  # top level
    B4 = Block4(obj_name="block4", sub_analyses=[B1])  # top level

    # B2._connect()
    # ic(B2.connects)

    # B1.declare_design_vars(variables=["x"])
    # B4.test_combined_adjoint()

    # exit()

    B1.set_var_values(variables={"x": 0.75})

    # Construct the system
    sys = System(
        sys_name="sys",
        top_level_analysis_list=[B3, B4],
        log_name="flume.log",
        log_prefix="examples/simple/test_paropt",
    )

    # sys.declare_design_vars(global_var_name={"block1.x": {"ub": 1.0}})
    sys.declare_design_vars(global_var_name={"block1.x": {"lb": 0.5, "ub": 1.0}})
    sys.declare_objective(global_obj_name="block3.alpha", obj_scale=10.0)
    sys.declare_constraints(
        global_con_name={
            "block4.beta": {"rhs": 3.0, "direction": "leq"},
            "block1.y": {"direction": "both", "rhs": 1.25},
        }
    )
    sys.declare_foi(global_foi_name=["block2.z"])

    # Create the FlumeParOptInterface
    interface = FlumeParOptInterface(flume_sys=sys)

    # Construct the paropt problem for the Flume system
    paroptprob = interface.construct_paropt_problem()
    options = interface.get_paropt_default_options(
        output_prefix="examples/simple/test_paropt"
    )

    # # Graph the system
    # graph = sys.graph_network()

    # graph.view()

    # exit()

    # Check the gradients for the system
    for i in range(5):
        paroptprob.checkGradients(1e-6)

    # Perform the optimization with ParOpt
    opt = ParOpt.Optimizer(paroptprob, options)
    opt.optimize()

    # Extract the optimized point
    x, z, zw, zl, zu = opt.getOptimizedPoint()

    ic(np.array(x))

    # interface.getVarsAndBounds(x, lb, ub)
    # ic(interface.ndvs)

    exit()

    # Check the full analysis list
    ic(sys.full_analysis_list)

    # Test the construction of the objective, constraint, design var, and FOI info
    sys.declare_objective(global_obj_name="block4.alpha")
    print("Objective Info:")
    ic(sys.obj_analysis, sys.obj_local_name)

    sys.declare_constraints(global_con_name={"block2.z": 0.0, "block4.beta": 1.0})
    print("\nConstraint Info:")
    for key in sys.con_info:
        ic(key, sys.con_info[key])

    sys.declare_design_vars(global_var_name={"block1.x": {"lb": 0.0, "ub": 1.0}})
    print("\nDesign Variable Info:")
    for key in sys.design_vars_info:
        ic(key, sys.design_vars_info[key])

    sys.declare_foi(global_foi_name=["block1.q"])
    print("\nFunctions of Interest Info:")
    for key in sys.foi:
        ic(key, sys.foi[key])

    exit()

    # Graph the system
    graph = sys.graph_network()

    graph.view()

    # plt.subplot(111)

    # nx.draw(graph, with_labels=True, font_weight="bold")

    # plt.show()
