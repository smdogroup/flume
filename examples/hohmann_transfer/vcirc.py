# from base_classes.analysis_base import AnalysisBase
from base_classes.analysis_base import AnalysisBase
from base_classes.state import State
import numpy as np


class VCircAnalysis(AnalysisBase):
    """
    Analysis class that computes the velcity associated with the circular orbit defined by the input radius and graviational parameter.
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
        r_var = State(value=1.0, desc="Radius from central body (km)", source=self)

        mu_var = State(
            value=1.0,
            desc="Gravitational parameter of central body (km**3/s**2)",
            source=self,
        )

        # Construct the variables dictionary
        self.variables = {"r": r_var, "mu": mu_var}

    def _analyze(self):
        """
        Private analysis method for the object that computes all of the outputs with the provided input variables and parameters to the object.
        """

        # Compute the circular velocity for the orbit
        v_c = np.sqrt(self.variables["mu"].value / self.variables["r"].value)

        # Update the attribute to reflect that the object has been analyzed
        self.analyzed = True

        # Assign the outputs to the outputs dictionary
        self.outputs = {}

        self.outputs["v_c"] = State(
            value=v_c,
            desc="Circular orbit velocity for the given radius and gravitational parameter (km/s)",
            source=self,
        )

        return

    def _analyze_adjoint(self):
        """
        Private adjoint analysis method for the object that computes the derivative values for the input variables via the adjoint method.
        """

        # Extract the derivative values for the outputs
        v_cb = self.outputs["v_c"].deriv

        # Compute the derivative contribution to r from v_c
        rb = (
            v_cb
            * 0.5
            * (self.variables["mu"].value / self.variables["r"].value) ** (-1 / 2)
            * -self.variables["mu"].value
            / self.variables["r"].value ** 2
        )

        # Compute the derivative contribution to mu from v_c
        mub = (
            v_cb
            * 0.5
            * (self.variables["mu"].value / self.variables["r"].value) ** (-1 / 2)
            / self.variables["r"].value
        )

        # Set the derivative values for the variables
        self.variables["r"].set_deriv_value(deriv_val=rb)
        self.variables["mu"].set_deriv_value(deriv_val=mub)

        # Update the attribute for the adjoint analysis
        self.adjoint_analyzed = True

        return


if __name__ == "__main__":
    # Create a VCirc analysis object
    vcirc = VCircAnalysis(obj_name="test")

    # Declare design variables for the object
    vcirc.declare_design_vars(variables=["mu", "r"])

    # Test the isolated adjoint implementation
    vcirc.test_isolated_adjoint(debug_print=True)
