from flume.base_classes.analysis import Analysis
from flume.base_classes.state import State
import numpy as np
from typing import List, Union
from .delta_v import DeltaVAnalysis


class TotalDeltaV(Analysis):
    def __init__(
        self,
        obj_name: str,
        sub_analyses=List[Union[DeltaVAnalysis, DeltaVAnalysis]],
        **kwargs,
    ):
        """
        DOCS:
        """

        # Set the default parameters
        self.default_parameters = {}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the default states for the variables
        dv1_var = State(value=0.0, desc="Delta V for the first maneuver", source=self)
        dv2_var = State(value=0.0, desc="Delta V for the second maneuver", source=self)

        self.variables = {"delta_v1": dv1_var, "delta_v2": dv2_var}

        return

    def _analyze(self):
        """
        Compute the total delta V from the delta V values for the two maneuvers.
        """

        # Extract the variable values
        dv1 = self.variables["delta_v1"].value
        dv2 = self.variables["delta_v2"].value

        # Compute the total delta_v
        dv_tot = dv1 + dv2

        # Update the attribute to reflect that the object has been analyzed
        self.analyzed = True

        # Store the outputs in the output dictionary
        self.outputs = {}

        self.outputs["delta_v_tot"] = State(
            value=dv_tot, desc="Total delta V for the Hohmann transfer", source=self
        )

        return

    def _analyze_adjoint(self):
        """
        Propagate the derivatives of the total delta V back to the individual delta V's.
        """

        # Extract the output/variable derivatives
        dv_totb = self.outputs["delta_v_tot"].deriv

        dv1b = self.variables["delta_v1"].deriv
        dv2b = self.variables["delta_v2"].deriv

        # Update the derivative values for the variables
        dv1b += dv_totb
        dv2b += dv_totb

        # Update the attribute for the adjoint analysis
        self.adjoint_analyzed = True

        # Set the derivative values
        self.variables["delta_v1"].set_deriv_value(deriv_val=dv1b)
        self.variables["delta_v2"].set_deriv_value(deriv_val=dv2b)

        return
