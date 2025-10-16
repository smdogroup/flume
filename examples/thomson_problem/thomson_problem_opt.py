from flume.base_classes.system import System
from icecream import ic
import numpy as np
from examples.thomson_problem.thomson_problem_classes import (
    ParticleConstraints,
    ParticlePositions,
    PotentialEnergy,
)
from math import pi
import argparse
from flume.interfaces.scipy_interface import FlumeScipyInterface
import matplotlib.pyplot as plt

if __name__ == "__main__":

    # Define an argument paraser
    parser = argparse.ArgumentParser(
        description="Program that performs the optimization of the Thomson problem, which seeks to find the minimum electrostatic potential energy configuration for n_p electrons constrained to the surface of a sphere."
    )

    # Add an argument that specifies the number of particles
    parser.add_argument(
        "--num-particles",
        "-np",
        type=int,
        default=3,
        help="Number of particles to consider for the system",
    )

    # Add an argument that defines the number of iterations
    parser.add_argument(
        "--max-iterations",
        "-maxit",
        type=int,
        default=100,
        help="Maximum number of iterations for the optimization problem",
    )

    # Add an argument that specifies whether the particle positions should be plotted
    parser.add_argument(
        "--plot-particles",
        "-plot",
        action="store_true",
        help="Argument that, when provided, triggers the script to plot the particle positions on the sphere after the optimization concludes.",
    )

    # Parse the arguments
    args = parser.parse_args()

    n_p = args.num_particles

    # Construct the analysis objects for the system
    positions = ParticlePositions(obj_name="positions", sub_analyses=[], n_p=n_p)

    energy = PotentialEnergy(obj_name="energy", sub_analyses=[positions], n_p=n_p)

    cons = ParticleConstraints(obj_name="cons", sub_analyses=[positions], n_p=n_p)

    # Construct the system
    sys = System(
        sys_name="thomson_problem",
        top_level_analysis_list=[energy, cons],
        log_name=f"flume_{n_p}.log",
        log_prefix="examples/thomson_problem",
    )

    # Declare the design variables for the system
    sys.declare_design_vars(
        global_var_name={
            "positions.theta": {"lb": -pi, "ub": pi},
            "positions.phi": {"lb": 0.0, "ub": 2 * pi},
        }
    )

    # Declare the objective
    sys.declare_objective(global_obj_name="energy.f")

    # Declare the constraints
    sys.declare_constraints(
        global_con_name={"cons.c": {"direction": "both", "rhs": 0.0}}
    )

    # Construct the Scipy interface
    interface = FlumeScipyInterface(flume_sys=sys)

    # Set random positions for x, y, z to start
    theta0 = np.random.uniform(size=n_p)
    phi0 = np.random.uniform(size=n_p)

    # Set the initial guess using theta0 and phi0
    init_global_vars = {"positions.theta": theta0, "positions.phi": phi0}
    initial_point = interface.set_initial_point(initial_global_vars=init_global_vars)

    # Optimize the problem with SciPy minimize (SLSQP does not work with this implementation because the subproblem is singular, so use COBYQA instead)
    maxit = max(100, n_p * 20)
    x, res = interface.optimize_system(
        x0=initial_point, method="trust-constr", maxit=maxit
    )

    # Display the result
    ic(res)

    # Check that the potential energy at the final point matches the expected value
    obj_val = res.fun

    def compute_rel_error(obj_val, obj_star):
        rel_error = abs(obj_star - obj_val) / obj_star

        print("\n%15s %15s" % ("Optimized PE:", f"{obj_val:.6f}"))
        print("%15s %15s" % ("Expected PE:", f"{obj_star:.6f}"))
        print("%15s %15s" % ("Rel. Error:", f"{rel_error:.6e}"))

        return None

    # Check the objective function value against the known solution if n_p has a known value
    if n_p == 2:
        obj_star = 0.5
        compute_rel_error(obj_val, obj_star)
    elif n_p == 3:
        obj_star = 1.732050808
        compute_rel_error(obj_val, obj_star)
    elif n_p == 12:
        obj_star = 49.165253058
        compute_rel_error(obj_val, obj_star)
    else:
        print(
            f"\nPotential energy for n_p = {n_p} = {obj_val:.4f}. No exact solution to compare against."
        )

    # Plot the particles on the sphere, if requested
    if args.plot_particles:
        # Create the figure
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")

        # Extract the values of particle positions
        outs = positions.get_output_values()
        x = outs["x"]
        y = outs["y"]
        z = outs["z"]

        # Plot the data for the sphere
        theta_sphere = np.linspace(-np.pi, np.pi, 20)
        phi_sphere = np.linspace(0, 2 * np.pi, 20)

        u, v = np.meshgrid(theta_sphere, phi_sphere)

        # Plot the surface of the sphere (unit radius)
        x_sphere = np.sin(u) * np.cos(v)
        y_sphere = np.sin(u) * np.sin(v)
        z_sphere = np.cos(u)

        ax.plot_surface(
            x_sphere,
            y_sphere,
            z_sphere,
            color="#f1535b",
            alpha=0.2,
            antialiased=True,
            shade=False,
            edgecolor="#989898",
            linewidth=0.5,
        )

        # Plot the locations of the points
        ax.plot(
            x,
            y,
            z,
            color="#0098cf",
            marker="o",
            markeredgecolor="None",
            linestyle="None",
            markersize=8,
        )

        # Set the axis labels
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.axis("equal")
        ax.grid(False)

        plt.show()
