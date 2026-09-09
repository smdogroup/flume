import openmdao.api as om
from flume.base_classes.om_wrapper import OpenMDAOGroupAnalysis
from components_and_group import BeamGroup
from flume.base_classes.system import System
from flume.interfaces.scipy_interface import FlumeScipyInterface
import numpy as np
from icecream import ic
import matplotlib.pyplot as plt

if __name__ == "__main__":

    # Set the values for the beam group/model
    E = 1.0
    L = 1.0
    b = 0.1
    volume = 0.01

    num_elements = 50

    # Initialize the OpenMDAO group for the beam model
    beam_group = BeamGroup(E=E, L=L, b=b, volume=volume, num_elements=num_elements)

    # Construct the Flume object that wraps the OpenMDAO group
    beam_flume = OpenMDAOGroupAnalysis(
        om_group=beam_group,
        group_variable_names=["beam_group.h"],
        group_output_names=[
            "beam_group.volume_comp.volume",
            "beam_group.compliance_comp.compliance",
        ],
        obj_name="beam_group",
    )

    # Construct the Flume System
    sys = System(sys_name="beam_opt", top_level_analysis_list=[beam_flume])

    # Declare the design variables
    sys.declare_design_vars(global_var_name={"beam_group.h": {"lb": 1e-2, "ub": 10.0}})

    # Declare the objective function for the beam
    obj_scale = 1e-5
    sys.declare_objective(
        global_obj_name="beam_group.compliance_comp.compliance", obj_scale=obj_scale
    )

    # Declare the constraint for the beam
    sys.declare_constraints(
        global_con_name={
            "beam_group.volume_comp.volume": {"direction": "both", "rhs": volume}
        }
    )

    # Construct the FlumeScipyInterface
    interface = FlumeScipyInterface(flume_sys=sys, callback=None)

    # Set the initial point for the optimization
    h0 = np.random.uniform(low=0.05, high=0.15, size=num_elements)
    x0 = interface.set_initial_point(initial_global_vars={"beam_group.h": h0})

    # Perform the optimization
    xstar, res = interface.optimize_system(
        x0=x0, options=None, method="SLSQP", maxit=300
    )

    ic(res)

    # Extract the optimal value of the compliance
    c = res.fun / obj_scale

    # Plot the optimized solution against the solution from OpenMDAO to compare
    cstar = 23762.153677294387

    rel_error = abs(cstar - c) / c

    print("\n%22s %15s" % ("Optimal Compliance:", f"{c:.6f}"))
    print("%22s %15s" % ("Expected Compliance:", f"{cstar:.6f}"))
    print("%22s %15s" % ("Rel. Error:", f"{rel_error:.6e}"))

    fig, ax = plt.subplots(1, 1)
    x = np.linspace(0.0, L)
    ax.plot(
        x,
        xstar,
        color="#0098cf",
        linewidth=1.2,
        marker="o",
        markeredgecolor="None",
        label="Flume",
    )

    hstar = np.array(
        [
            0.14915747,
            0.14764318,
            0.14611303,
            0.14456710,
            0.14300483,
            0.14142408,
            0.13982622,
            0.13820997,
            0.13657415,
            0.13491848,
            0.13324256,
            0.13154536,
            0.12982559,
            0.12808303,
            0.12631664,
            0.12452488,
            0.12270693,
            0.12086161,
            0.11898795,
            0.11708423,
            0.11514929,
            0.11318075,
            0.11117751,
            0.10913767,
            0.10705892,
            0.10493900,
            0.10277538,
            0.10056522,
            0.09830538,
            0.09599246,
            0.09362258,
            0.09119082,
            0.08869253,
            0.08612198,
            0.08347228,
            0.08073578,
            0.07790314,
            0.07496381,
            0.07190458,
            0.06870939,
            0.06535830,
            0.06182635,
            0.05808044,
            0.05407648,
            0.04975297,
            0.04501847,
            0.03972915,
            0.03363155,
            0.02620193,
            0.01610861,
        ]
    )

    diff_norm = (np.sum((xstar - hstar) ** 2)) ** 0.5
    ic(diff_norm)

    ax.plot(
        x,
        hstar,
        color="#f1535b",
        linewidth=1.2,
        marker="^",
        markeredgecolor="None",
        linestyle="None",
        label="OpenMDAO",
    )

    ax.set_title(f"Norm of Difference = {diff_norm:.4e}")
    ax.set_xlabel("X")
    ax.set_ylabel("Thickness Distribution (h)")
    ax.legend(loc="lower left")

    plt.show()
