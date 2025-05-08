from flume.base_classes.analysis import Analysis
from flume.base_classes.state import State
import numpy as np
from icecream import ic
from typing import List, Union
from math import exp


class FirstBlock(Analysis):

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

        # Compute y
        y = x**2

        # Update analyzed attribute
        self.analyzed = True

        # Store the outputs
        self.outputs = {}

        self.outputs["y"] = State(value=y, desc="y value", source=self)

    def _analyze_adjoint(self):
        # Extract derivative of the output
        yb = self.outputs["y"].deriv

        # Extract the derivative value for x
        xb = self.variables["x"].deriv

        # Compute xb from y
        x = self.variables["x"].value
        xb += 2 * x * yb

        # Update the analyzed adjoint attribute
        self.adjoint_analyzed = True

        # Update the derivative values for the varialbes
        self.variables["x"].set_deriv_value(deriv_val=xb)

        return


class ObjTerm1(Analysis):

    def __init__(self, obj_name: str, sub_analyses=List[FirstBlock], **kwargs):
        # Set default parameters
        self.default_parameters = {}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the defualt State for the variable
        y_var = State(value=1.0, desc="y value", source=self)

        # Construct vars dictionary
        self.variables = {"y": y_var}

        return

    def _analyze(self):
        # Extract variable values
        y = self.variables["y"].value

        # Compute f
        f = 2 * y

        # Update analyzed attribute
        self.analyzed = True

        # Store the outputs
        self.outputs = {}

        self.outputs["f"] = State(value=f, desc="f value", source=self)

        return

    def _analyze_adjoint(self):
        # Extract the output derivative
        fb = self.outputs["f"].deriv

        # Extract the variable derivative
        yb = self.variables["y"].deriv

        # Compute the contribution to yb from f
        yb += 2 * fb

        # Update the analyzed adjoint attribute
        self.adjoint_analyzed = True

        # Update the derivative values for the varialbes
        self.variables["y"].set_deriv_value(deriv_val=yb)

        return


class ObjTerm2(Analysis):

    def __init__(self, obj_name: str, sub_analyses=List[FirstBlock], **kwargs):
        # Set default parameters
        self.default_parameters = {}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the defualt State for the variable
        y_var = State(value=1.0, desc="y value", source=self)

        # Construct vars dictionary
        self.variables = {"y": y_var}

        return

    def _analyze(self):
        # Extract variable values
        y = self.variables["y"].value

        # Compute g
        g = np.exp(y).item()

        # Update analyzed attribute
        self.analyzed = True

        # Store the outputs
        self.outputs = {}

        self.outputs["g"] = State(value=g, desc="g value", source=self)

        return

    def _analyze_adjoint(self):
        # Extract the output derivative
        gb = self.outputs["g"].deriv

        # Extract the variable derivative
        yb = self.variables["y"].deriv

        # Compute the contribution to yb from f
        y = self.variables["y"].value
        yb += exp(y) * gb

        # Update the analyzed adjoint attribute
        self.adjoint_analyzed = True

        # Update the derivative values for the varialbes
        self.variables["y"].set_deriv_value(deriv_val=yb)

        return


class MultiObjective(Analysis):

    def __init__(
        self, obj_name: str, sub_analyses=List[Union[ObjTerm1, ObjTerm2]], **kwargs
    ):

        # Set default parameters
        self.default_parameters = {"w": 0.5}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the default State for the variables
        f_var = State(value=1.0, desc="f state value", source=self)

        g_var = State(value=1.0, desc="g state value", source=self)

        # Construct variables dictionary
        self.variables = {"f": f_var, "g": g_var}

        return

    def _analyze(self):
        # Extract variable values
        f = self.variables["f"].value
        g = self.variables["g"].value
        w = self.parameters["w"]

        # Compute the multi objective value
        # J = w * f + (1 - w) * g
        # J = w * f
        J = (1 - w) * g

        # Update the analyzed attribute
        self.analyzed = True

        # Store the outputs
        self.outputs = {}

        self.outputs["J"] = State(
            value=J, desc="Multi-objective function value", source=self
        )

        return

    def _analyze_adjoint(self):
        # Extract output derivative
        Jb = self.outputs["J"].deriv

        # Extract the variable derivatives
        fb = self.variables["f"].deriv
        gb = self.variables["g"].deriv

        # Extract the weight value
        w = self.parameters["w"]

        # Compute the contribution to fb from J
        fb += Jb * w

        # Compute the contribution to gb from J
        gb += Jb * (1 - w)

        # Update the analyzed adjoint attribute
        self.adjoint_analyzed = True

        # Update the derivative values for the variables
        # self.variables["f"].set_deriv_value(deriv_val=fb)
        self.variables["g"].set_deriv_value(deriv_val=gb)

        return


if __name__ == "__main__":

    # Construct the four analysis blocks
    first = FirstBlock(obj_name="first")

    term1 = ObjTerm1(obj_name="term1", sub_analyses=[first])
    term2 = ObjTerm2(obj_name="term2", sub_analyses=[first])

    w = 0.5
    multi = MultiObjective(obj_name="multi", sub_analyses=[term1, term2], w=w)

    # Test the isolated adjoints
    test_isolated_adjoints = False
    if test_isolated_adjoints:
        first.declare_design_vars(variables=["x"])
        first.test_isolated_adjoint(method="cs")

        term1.declare_design_vars(variables=["y"])
        term1.test_isolated_adjoint(method="cs")

        term2.declare_design_vars(variables=["y"])
        term2.test_isolated_adjoint(method="cs")

        multi.declare_design_vars(variables=["f", "g"])
        multi.test_isolated_adjoint(method="cs")

    # Test the isolated adjoing for the system
    first.declare_design_vars(variables=["x"])
    multi.test_combined_adjoint(method="cs", debug_print=False)
