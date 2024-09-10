from base_classes.analysis_base import Analysis
from base_classes.state import State
from examples.hohmann_transfer.vcirc import VCircAnalysis
from examples.hohmann_transfer.transfer_orbit import TransferOrbitAnalysis
from examples.hohmann_transfer.delta_v import DeltaVAnalysis
import numpy as np


class HohmannLEOtoGEO:
    """
    Class that combines the different analyses needed to perform an analysis of a Hohmann transfer from LEO go GEO.
    """

    def __init__(
        self,
        leo: VCircAnalysis,
        geo: VCircAnalysis,
        transfer: TransferOrbitAnalysis,
        dv1: DeltaVAnalysis,
        dv2: DeltaVAnalysis,
    ):
        """
        Class that combines analyses to perform Hohmann transfer analysis from LEO to GEO.

        Parameters
        ----------
        leo : VCircAnalysis
            Instance of VCircAnalysis that defines the LEO orbit
        geo : VCircAnalysis
            Instance of VCircAnalysis that defines the GEO orbit
        transfer : TransferOrbitAnalysis
            Instance of TransferOrbitAnalysis that defines the parameters for the transfer orbit between the two circular orbits.
        dv1 : DeltaVAanlysis
            Instance of DeltaVAnalysis that computes the first velocity change required based on velocity and inclination data
        dv2 : DeltaVAanlysis
            Instance of DeltaVAnalysis that computes the second velocity change required based on velocity and inclination data
        """

        # Store the different objects as attributes
        self.leo = leo
        self.geo = geo
        self.transfer = transfer
        self.dv1 = dv1
        self.dv2 = dv2

        return

    def perform_hohmann(self):
        """
        Performs the Hohmann transfer analysis with input objects.

        Returns
        -------
        delta_v_total : float
            Total velocity change (km/s) required for the Hohmann transfer (sum of the delta v's for the two maneuvers)
        dinc_tot : float
            Total inclination change (rad) for the Hohmann transfer (sum of the individual inclination changes for the two manuevers)
        """

        # Analyze the LEO and GEO circular orbits
        self.leo.analyze()
        self.geo.analyze()

        # Perform the transfer orbit analysis
        self.transfer.analyze()

        # Extract the outputs from each of the analyses
        leo_outs = self.leo.get_output_values()
        geo_outs = self.geo.get_output_values()
        transfer_outs = self.transfer.get_output_values()

        # Set the variable values for the first delta_v analysis (i.e. LEO to transfer)
        self.dv1.set_var_values(
            variables={"v1": leo_outs["v_c"], "v2": transfer_outs["vp"]}
        )

        # Set the variable values for the second delta_v analysis (i.e. transfer to GEO)
        self.dv2.set_var_values(
            variables={"v1": transfer_outs["va"], "v2": geo_outs["v_c"]}
        )

        # Analyze the two delta_v analyses
        self.dv1.analyze()
        self.dv2.analyze()

        # Extract the delta_v values from the two analyses
        dv1_outs = self.dv1.get_output_values()
        delta_v1 = dv1_outs["delta_v"]

        dv2_outs = self.dv2.get_output_values()
        delta_v2 = dv2_outs["delta_v"]

        # Compute the total delta_v
        delta_v_total = delta_v1 + delta_v2

        # Compute the total inclination change
        dinc_tot = self.dv1.variables["dinc"].value + self.dv2.variables["dinc"].value

        return delta_v_total, dinc_tot


if __name__ == "__main__":
    # Define the inputs for the system
    mu = 398600.4418
    r1 = 6778.0
    r2 = 42164.0
    dinc1 = 0.0
    dinc2 = 28.5 * np.pi / 180

    # Create the analysis objects for the various parts of the Hohmann transfer maneuver
    leo = VCircAnalysis(obj_name="leo")
    leo.set_var_values(variables={"mu": mu, "r": r1})

    geo = VCircAnalysis(obj_name="geo")
    geo.set_var_values(variables={"mu": mu, "r": r2})

    transfer = TransferOrbitAnalysis(obj_name="transfer", mu=mu)
    transfer.set_var_values(variables={"ra": r2, "rp": r1})

    dv1 = DeltaVAnalysis(obj_name="dv1")
    dv1.set_var_values(variables={"dinc": dinc1})

    dv2 = DeltaVAnalysis(obj_name="dv2")
    dv2.set_var_values(variables={"dinc": dinc2})

    # Construct the HohmannLEOtoGEO system
    hohmann = HohmannLEOtoGEO(leo=leo, geo=geo, transfer=transfer, dv1=dv1, dv2=dv2)

    # Perform the Hohmann transfer analysis
    dv_tot, dinc_tot = hohmann.perform_hohmann()

    print("Delta-V (km/s):", dv_tot)
    print("Inclination change split (deg):", dinc_tot * 180 / np.pi)
