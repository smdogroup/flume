from base_classes.analysis import Analysis
from base_classes.state import State
import numpy as np


class DeltaVAnalysis(Analysis):
    """
    Analysis class that computes the change in velocity required to perform a maneuver given the magnitudes of the velocities before and after the manuever and the inclination change between the orbits.
    """

    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):
        """
        Analysis class for computing change in velocity of two orbits.

        Parameters
        ----------
        obj_name : str
            Name for the analysis object
        sub_analyses : list
            A list of any sub-analyses for the DeltaV analysis (nominally, this is empty)

        Keyword Arguments
        -----------------
        No input parameters for this object, so no **kwargs.

        """

        # Set the default parameters for the object
        self.default_parameters = {}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the default States for the variables
        v1_var = State(value=1.0, desc="Initial velocity (km/s)", source=self)

        v2_var = State(value=1.0, desc="Final velocity (km/s)", source=self)

        dinc_var = State(value=1.0, desc="Plane change (rad)", source=self)

        # Construct the variables dictionary
        self.variables = {"v1": v1_var, "v2": v2_var, "dinc": dinc_var}

        return

    def _analyze(self):
        """
        Private analysis method for the object that computes all of the outputs with the provided input variables and parameters to the object.
        """

        # Compute the change in velocity
        delta_v = np.sqrt(
            self.variables["v1"].value ** 2
            + self.variables["v2"].value ** 2
            - 2
            * self.variables["v1"].value
            * self.variables["v2"].value
            * np.cos(self.variables["dinc"].value)
        )

        # Update the attribute to reflect that the object has been analyzed
        self.analyzed = True

        # Store the output in the output dictionary
        self.outputs = {}

        self.outputs["delta_v"] = State(
            value=delta_v, desc="Change in velocity (km/s)", source=self
        )

        return

    def _analyze_adjoint(self):
        """
        Private adjoint analysis method for the object that computes the derivative values for the input variables via the adjoint method.
        """

        # Extract the existing derivative values for the outputs and variables
        delta_vb = self.outputs["delta_v"].deriv

        # Compute the derivative contribution to v1 from delta_v
        v1b = (
            delta_vb
            * 0.5
            / self.outputs["delta_v"].value
            * (
                2 * self.variables["v1"].value
                - 2 * self.variables["v2"].value * np.cos(self.variables["dinc"].value)
            )
        )

        v2b = (
            delta_vb
            * 0.5
            / self.outputs["delta_v"].value
            * (
                2 * self.variables["v2"].value
                - 2 * self.variables["v1"].value * np.cos(self.variables["dinc"].value)
            )
        )

        dincb = (
            delta_vb
            * 0.5
            / self.outputs["delta_v"].value
            * (
                2
                * self.variables["v1"].value
                * self.variables["v2"].value
                * np.sin(self.variables["dinc"].value)
            )
        )

        # Set the derivative values for the variables
        self.variables["v1"].set_deriv_value(deriv_val=v1b)
        self.variables["v2"].set_deriv_value(deriv_val=v2b)
        self.variables["dinc"].set_deriv_value(deriv_val=dincb)

        # Update the attribute for the adjoint analysis
        self.adjoint_analyzed = True

        return


if __name__ == "__main__":
    # Create a deltav analysis object
    deltav = DeltaVAnalysis(obj_name="test")

    # Declare design variables for the object
    deltav.declare_design_vars(variables=["v1", "v2", "dinc"])

    # Test the isolated adjoint implementation
    deltav.test_isolated_adjoint(debug_print=True)
