# Description: the example outlined in this file is adapted from the OpenMDAO documentation "Hohmann Transfer Example - Optimizing a Spacecraft Manuever" (https://openmdao.org/newdocs/versions/latest/examples/hohmann_transfer/hohmann_transfer.html). This file defines the run script that combines the individual Flume Analysis classes, located within the other files in this directory, into a System. The FlumeScipyInterface is utilized to perform the optimization, as well.

from flume.base_classes.system import System
from examples.hohmann_transfer.vcirc import VCircAnalysis
from examples.hohmann_transfer.transfer_orbit import TransferOrbitAnalysis
from examples.hohmann_transfer.delta_v import DeltaVAnalysis
import numpy as np
from .total_delta_v import TotalDeltaV
from .inclination_change import TotalInclinationChange
from .independents import Independents
from flume.interfaces.scipy_interface import FlumeScipyInterface
from icecream import ic


if __name__ == "__main__":

    # Define the inputs for the system
    mu = 398600.4418
    r1 = 6778.0
    r2 = 42164.0
    dinc1 = 0.0
    dinc2 = 28.5 * np.pi / 180

    # Construct the various analysis objects
    indeps = Independents(obj_name="indeps")
    indeps.set_var_values(
        variables={
            "mu_dv": mu,
            "r1_dv": r1,
            "r2_dv": r2,
            "dinc1_dv": dinc1,
            "dinc2_dv": dinc2,
        }
    )

    leo = VCircAnalysis(
        obj_name="leo",
        sub_analyses=[indeps],
        variable_aliases={"r": "r1", "mu": "mu"},
        output_aliases={"v_c": "v_cl"},
    )

    geo = VCircAnalysis(
        obj_name="geo",
        sub_analyses=[indeps],
        variable_aliases={"r": "r2", "mu": "mu"},
        output_aliases={"v_c": "v_cg"},
    )

    transfer = TransferOrbitAnalysis(
        obj_name="transfer",
        sub_analyses=[indeps],
        variable_aliases={"ra": "r2", "rp": "r1", "mu": "mu"},
    )

    dv1 = DeltaVAnalysis(
        obj_name="dv1",
        sub_analyses=[indeps, leo, transfer],
        variable_aliases={"v1": "v_cl", "v2": "vp", "dinc": "dinc1"},
        output_aliases={"delta_v": "delta_v1"},
    )

    dv2 = DeltaVAnalysis(
        obj_name="dv2",
        sub_analyses=[indeps, geo, transfer],
        variable_aliases={"v1": "va", "v2": "v_cg", "dinc": "dinc2"},
        output_aliases={"delta_v": "delta_v2"},
    )

    obj = TotalDeltaV(obj_name="dv_tot", sub_analyses=[dv1, dv2])

    # To test the combined adjoint for the objective function, uncomment the following lines
    # indeps.declare_design_vars(variables=["dinc1_dv", "dinc2_dv"])
    # obj.test_combined_adjoint(method="cs")

    con = TotalInclinationChange(obj_name="dinc_tot", sub_analyses=[indeps])

    # To test the combined adjoint for the constraint function, uncomment the following lines
    # indeps.declare_design_vars(variables=["dinc1_dv", "dinc2_dv"])
    # con.test_combined_adjoint(method="cs")

    # Setup the Flume system
    sys = System(
        sys_name="Hohmann_Transfer",
        top_level_analysis_list=[obj, con],
        log_name="flume.log",
        log_prefix="examples/hohmann_transfer",
    )

    # Graph the network defined by the System and save to the output directory, which is used to visualize the structure of the System
    graph = sys.graph_network(
        filename="HohmannSystem", output_directory="examples/hohmann_transfer"
    )

    # Declare the objective function value
    sys.declare_objective(global_obj_name="dv_tot.delta_v_tot")

    # Declare the equality constraint for the total inclination change
    sys.declare_constraints(
        global_con_name={
            "dinc_tot.dinc_tot": {"direction": "both", "rhs": 28.5 * np.pi / 180}
        }
    )

    # Declare the design variables for the System
    sys.declare_design_vars(
        global_var_name={
            "indeps.dinc1_dv": {"lower": 0.0, "upper": 28.5 * np.pi / 180},
            "indeps.dinc2_dv": {"lower": 0.0, "upper": 28.5 * np.pi / 180},
        }
    )

    # Construct the FlumeScipyInterface
    interface = FlumeScipyInterface(flume_sys=sys)

    # Set the initial point for the optimizer
    x0 = interface.set_initial_point(
        initial_global_vars={
            "indeps.dinc1_dv": 0.0 * np.pi / 180,
            "indeps.dinc2_dv": 28.5 * np.pi / 180,
        }
    )

    # Optimize the system and print the results
    xstar, res = interface.optimize_system(x0=x0, maxit=100)

    ic(res)
    ic(xstar * 180 / np.pi)

    obj_outs = obj.get_output_values(outputs=["delta_v_tot"])
    dv_tot = obj_outs["delta_v_tot"]
    con_outs = con.get_output_values(outputs=["dinc_tot"])
    dinc_tot = con_outs["dinc_tot"]

    print("\nTotal Delta-V (km/s):", dv_tot)
    print("Inclination change split (deg):", dinc_tot * 180 / np.pi)
