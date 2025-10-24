# Description: the example outlined in this file is adapted from the OpenMDAO documentation "Optimizing the Thickness Distribution of a Cantilever Beam Using the Adjoint Method" (https://openmdao.org/newdocs/versions/latest/examples/beam_optimization_example.html). This file defines the run script, which combines the individual Flume Analysis classes defined within the directory for this example. Optimization is performed using the FlumeScipyInterface, and a graphical comparison is made with the solution from OpenMDAO.

from flume.base_classes.system import System
from flume.interfaces.scipy_interface import FlumeScipyInterface
from .independents import Independents
from .compliance import Compliance
from .inertia import MomentofInertia
from .linear_solve import LinearStaticSolve
from .volume import BeamVolume
import numpy as np
from icecream import ic
import matplotlib.pyplot as plt


if __name__ == "__main__":

    # Set the parameter values
    L = 1.0
    b = 0.1
    E = 1.0

    volume_constraint = 0.01
    nelems = 50

    # Compute the number of nodes
    nnodes = nelems + 1

    # Construct the force vector (point load in the vertical direction at the final node)
    force_vector = np.zeros(2 * nnodes)
    force_vector[-2] = -1.0

    # Construct the Independents object
    indeps = Independents(obj_name="indeps", sub_analyses=[], nelems=nelems)

    # Construct the MomentofInertia object
    inertia = MomentofInertia(
        obj_name="inertia", sub_analyses=[indeps], b=b, nelems=nelems
    )

    # Construct the LinearStaticSolve object
    linear_solve = LinearStaticSolve(
        obj_name="linear_solve",
        sub_analyses=[inertia],
        nelems=nelems,
        E=E,
        L=L,
        f=force_vector,
    )

    # Construct the Compliance object
    compliance = Compliance(
        obj_name="compliance",
        sub_analyses=[linear_solve],
        nelems=nelems,
        f=force_vector,
    )

    # Construct the BeamVolume object
    volume = BeamVolume(
        obj_name="volume", sub_analyses=[indeps], nelems=nelems, b=b, L=L
    )

    # Construct the System for the optimization problem
    sys = System(
        sys_name="BeamThicknessOptimization",
        top_level_analysis_list=[compliance, volume],
        log_name="flume.log",
        log_prefix="examples/cantilever_beam",
    )

    # Declare the objective function for the beam
    obj_scale = 1e5
    sys.declare_objective(global_obj_name="compliance.c", obj_scale=1e-5)

    # Declare the constraint for the beam
    sys.declare_constraints(
        global_con_name={"volume.V": {"direction": "both", "rhs": volume_constraint}}
    )

    # Declare the design variables
    sys.declare_design_vars(global_var_name={"indeps.h_dv": {"lb": 1e-2, "ub": 10.0}})

    # Graph the network defined by the System and save to the output directory, which is used to visualize the structure of the System
    graph = sys.graph_network(
        filename="CantileverBeam", output_directory="examples/cantilever_beam"
    )

    # Construct the FlumeScipyInterface
    interface = FlumeScipyInterface(flume_sys=sys, callback=None)

    # Set the initial point for the optimization
    h0 = np.random.uniform(low=0.05, high=0.15, size=nelems)
    x0 = interface.set_initial_point(initial_global_vars={"indeps.h_dv": h0})

    # Perform the optimization
    xstar, res = interface.optimize_system(
        x0=x0, options=None, method="SLSQP", maxit=300
    )

    ic(res)

    # Extract the optimal value of the compliance
    c = compliance.get_output_values()["c"]

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
