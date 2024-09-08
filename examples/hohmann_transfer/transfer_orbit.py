from base_classes.analysis_base import AnalysisBase
from base_classes.state import State
import numpy as np


class TransferOrbitAnalysis(AnalysisBase):
    """
    Analysis class that computes the velocity magnitudes at periapsis and apoapsis of an orbit for a given transfer orbit maneuver.
    """

    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):
        """
        Analysis clas for computing the velocities of the transfer orbit.

        Parameters
        ----------
        obj_name : str
            Name for the analysis object
        sub_analyses : list
            A list of any sub-analyses for the DeltaV analysis (nominally, this is empty)

        Keyword Arguments
        -----------------
        mu : float
            Graviational parameter of the central body for the orbit (km**3/s**2)
        """

        # Set the default parameters for the object
        self.default_parameters = {
            "mu": 398600.4418  # Graviational parameter of the central body (km**3/s**2)
        }

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the default States for the variables
        ra_var = State(value=42164.0, desc="Apoapsis radius (km)", source=self)

        rp_var = State(value=7000.0, desc="Periapsis radius (km)", source=self)

        # Construct the variables dictionary
        self.variables = {"ra": ra_var, "rp": rp_var}

        return

    def _analyze(self):
        """
        Private analysis method for the object that computes all of the outputs with the provided input variables and parameters to the object.
        """

        # Compute the semi-major axis of the transfer orbit
        self.a = (self.variables["ra"].value + self.variables["rp"].value) / 2.0

        # Compute the eccentricity of the transfer orbit
        self.e = (self.a - self.variables["rp"].value) / self.a

        # Compute the semilatus rectrum of the orbit
        self.p = self.a * (1.0 - self.e**2)

        # Compute the specific angular momentum of the orbit
        self.h = np.sqrt(self.p * self.parameters["mu"])

        # Compute the velocities at periapsis and apoapsis
        va = self.h / self.variables["ra"].value

        vp = self.h / self.variables["rp"].value

        # Update the attribute to reflect that the object has been analyzed
        self.analyzed = True

        # Store the outputs in the output dictionary
        self.outputs = {}

        self.outputs["va"] = State(
            value=va, desc="Velocity at apoapsis (km/s)", source=self
        )

        self.outputs["vp"] = State(
            value=vp, desc="Velocity at periapsis (km/s)", source=self
        )

        return

    def _analyze_adjoint(self):
        """
        Private adjoint analysis method for the object that computes the derivative values for the input variables via the adjoint method.
        """

        # Extract the exiting derivative values for the outputs and variables
        vab = self.outputs["va"].deriv
        vpb = self.outputs["vp"].deriv

        rab = self.variables["ra"].deriv
        rpb = self.variables["rp"].deriv

        # Compute the derivative contribution to self.h from va and vp
        hb = vab / self.variables["ra"].value
        hb += vpb / self.variables["rp"].value

        # Compute the contributions to ra from va
        rab += vab * -self.h / self.variables["ra"].value ** 2

        # Compute the contributions to rp from vp
        rpb += vpb * -self.h / self.variables["rp"].value ** 2

        # Compute the contributions to mu from h
        mub = hb * 0.5 * (self.parameters["mu"] * self.p) ** (-1 / 2) * self.p

        # Compute the contributions to p from h
        pb = (
            hb
            * 0.5
            * (self.parameters["mu"] * self.p) ** (-1 / 2)
            * self.parameters["mu"]
        )

        # Compute the contribution to a from p
        ab = pb * (1 - self.e**2)

        # Compute the contribution to e from p
        eb = pb * -2 * self.e * self.a

        # Compute the contribution to a from e
        ab += eb * self.variables["rp"].value / self.a**2

        # Compute the contribution to rp from e
        rpb += eb * -1 / self.a

        # Compute the contribution to ra from a
        rab += ab * 0.5

        # Compute the contribution to rp from a
        rpb += ab * 0.5

        # Set the derivative values for the variables
        self.variables["ra"].set_deriv_value(deriv_val=rab)
        self.variables["rp"].set_deriv_value(deriv_val=rpb)

        # Update the attribute for the adjoint analysis
        self.adjoint_analyzed = True

        return


if __name__ == "__main__":
    # Create the transfer orbit analysis object
    transfer = TransferOrbitAnalysis(obj_name="test")

    # Declare design variables for the object
    transfer.declare_design_vars(variables=["ra", "rp"])

    # Test the isolated adjoint implementation
    transfer.test_isolated_adjoint(debug_print=True)
