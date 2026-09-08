from paraboloid_component import Paraboloid
from flume.base_classes.state import State
from flume.base_classes.analysis import Analysis
from flume.base_classes.system import System
from flume.interfaces.scipy_interface import FlumeScipyInterface
import openmdao.api as om
from icecream import ic
from flume.base_classes.om_wrapper import OpenMDAOComponentAnalysis
import numpy as np


class DesignVars(Analysis):

    def __init__(self, obj_name: str, sub_analyses: list = [], **kwargs):

        # Declare no parameters
        self.default_parameters = {}

        # Super init
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Setup the variables

        x_var = State(value=0.0, desc="x value", source=self)
        y_var = State(value=0.0, desc="y value", source=self)

        self.variables = {"x_dv": x_var, "y_dv": y_var}

        return

    def _analyze(self):
        """
        Map dvs to output x and y.
        """

        # Extract design variable values and map to outputs
        x = self.variables["x_dv"].value
        y = self.variables["y_dv"].value

        # Assign to outputs
        self.outputs = {}

        self.outputs["x"] = State(value=x, desc="x value", source=self)
        self.outputs["y"] = State(value=y, desc="y value", source=self)

        return

    def _analyze_adjoint(self):
        """
        Map output derivs back to design variables
        """

        # Get the output deriv values and map to design variable derivs
        xb = self.outputs["x"].deriv
        yb = self.outputs["y"].deriv

        x_dvb = self.variables["x_dv"].deriv
        y_dvb = self.variables["y_dv"].deriv

        self.variables["x_dv"].set_deriv_value(xb + x_dvb)
        self.variables["y_dv"].set_deriv_value(yb + y_dvb)

        return


if __name__ == "__main__":

    # Construct the design variables object
    dv_obj = DesignVars(obj_name="dvs")

    # Construct the OpenMDAO Paraboloid component
    om_parab = Paraboloid()

    # Construct the Flume object that maps to the OpenMDAO Paraboloid component
    flume_parab = OpenMDAOComponentAnalysis(
        om_component=om_parab, obj_name="parab", sub_analyses=[dv_obj]
    )

    # Construct the System
    sys = System(
        sys_name="paraboloid_opt",
        top_level_analysis_list=[flume_parab],
    )

    # Declare the design variables
    sys.declare_design_vars(global_var_name={"dvs.x_dv": {}, "dvs.y_dv": {}})

    # Declare the objective
    sys.declare_objective(global_obj_name="parab.f_xy")

    # Setup the optimizer interface
    interface = FlumeScipyInterface(flume_sys=sys, callback=None)

    x0 = interface.set_initial_point(
        initial_global_vars={"dvs.x_dv": 3.0, "dvs.y_dv": -4.0}
    )

    xstar, res = interface.optimize_system(x0)

    ic(res)
    ic(xstar)
