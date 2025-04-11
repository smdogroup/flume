from paropt import ParOpt
from mpi4py import MPI
import os
from flume.base_classes.system import System
import numpy as np
from icecream import ic


class ParOptProb(ParOpt.Problem):
    def __init__(self, comm, prob) -> None:
        self.prob = prob
        super(ParOptProb, self).__init__(comm, nvars=prob.ndvs, ncon=prob.ncon)
        return

    def getVarsAndBounds(self, x, lb, ub):
        return self.prob.getVarsAndBounds(x, lb, ub)

    def evalObjCon(self, x):
        return self.prob.evalObjCon(x)

    def evalObjConGradient(self, x, g, A):
        return self.prob.evalObjConGradient(x, g, A)


class FlumeParOptInterface:

    def __init__(self, flume_sys: System, callback=None, update=None):
        """
        Creates an interface that is used to link an instance of a Flume System to a ParOpt optimization problem.

        Parameters
        ----------
        flume_sys : System
            Instance of a Flume System that represents the problem to be solved with ParOpt
        callback : callable function, default None
            This is a callable function that gets executed during every iteration of the evalObjCon method during optimization. See below for the structure of the function
        update : callable function, default None
            This is a callable function that gets executed at the start of every iteration of the evalObjCon method during optimization. Nominally, this is used to do parameter updates, such as for a continuation strategy

        Callback Function
        -----------------
        The callback function is an arbitrary function that the user writes, and it will be called during every evalObjCon execution. It is required that the function is setup such that it takes in two arguments:

        def user_callback(x, it_number):
            ...

        Here, the parameters are:
        x : current design variable info for the problem
        it_number : current iteration number for the system
        """

        # Store the flume system as an attribute
        self.flume_sys = flume_sys

        # Store the callback function
        self.callback = callback

        # Store the update function
        self.update = update

        # Store the string that declares the optimizer
        self.optimizer = "paropt"

        # Initialize the iteration counter
        self.it_counter = 0

        # Set the number of constraints
        if hasattr(self.flume_sys, "con_info"):
            self.ncon = len(self.flume_sys.con_info.keys())
        else:
            self.ncon = 0

        # Set the number of design variables and the indices for the ParOpt PVec x
        self._set_ndvs()

        return

    def _set_ndvs(self):
        """
        Gets the number of design variables for the system of interest.
        """

        ndvs = 0
        self.indices = {}
        # Loop through each of the keys in the design_vars_info dictionary
        for var in self.flume_sys.design_vars_info:

            # Extract the local name for the current variable
            local_name = self.flume_sys.design_vars_info[var]["local_name"]

            # Extract variable value for the current variable
            var_i = self.flume_sys.design_vars_info[var]["instance"].get_var_values(
                variables=[local_name]
            )[local_name]

            # Extract the number of design variables associated with the current variable
            if isinstance(var_i, float):
                ndvs_i = 1
            elif isinstance(var_i, np.ndarray):
                ndvs_i = np.size(var_i)

            # Store the start and end indices for the current design variable, which is used to set/extract the variables from the x PVec
            self.indices[var] = {"start": ndvs, "end": ndvs + ndvs_i}

            # Add the number of dvs to the total number of design variables
            ndvs += ndvs_i

        # Set the attribute for the number of dvs
        self.ndvs = ndvs

        return

    def getVarsAndBounds(self, x, lb, ub):
        """
        Get the variable values and set the bounds for the optimization problem.
        """

        # Check to make sure that the design variable info has been set, otherwise raise an error
        if not hasattr(self.flume_sys, "design_vars_info"):
            raise RuntimeError(
                f"The design variables have not yet been declared for the system named '{self.flume_sys.sys_name},' so the bounds cannot be set. Ensure that the function 'declare_design_vars' has been called."
            )

        # # Get the number of design variables for the problem, if necessary
        # if not hasattr(self, "ndvs"):
        #     self._set_ndvs()

        # Extract the design variable and bound values for each variable in the system
        tracked_dvs = 0
        for var in self.flume_sys.design_vars_info:
            # Extract the local name for the current variable
            local_name = self.flume_sys.design_vars_info[var]["local_name"]

            # Extract variable value for the current variable
            var_i = self.flume_sys.design_vars_info[var]["instance"].get_var_values(
                variables=[local_name]
            )[local_name]

            # ic(var_i)

            # # Extract the number of design variables associated with the current variable
            # if isinstance(var_i, float):
            #     ndvs_i = 1
            # elif isinstance(var_i, np.ndarray):
            #     ndvs_i = np.size(var_i)

            # Set the value in the x PVec
            start = self.indices[var]["start"]
            end = self.indices[var]["end"]
            x[start:end] = var_i

            # Set the lower bound value for the current variable, if it exists
            if "lb" in self.flume_sys.design_vars_info[var]:
                lb[start:end] = self.flume_sys.design_vars_info[var]["lb"]
            else:
                lb[start:end] = -1e30

            # Set the upper bound value for the current variable
            if "ub" in self.flume_sys.design_vars_info[var]:
                ub[start:end] = self.flume_sys.design_vars_info[var]["ub"]
            else:
                ub[start:end] = 1e30

        return

    def set_system_variables(self, x, it_counter):
        """
        Sets the variable values for the analysis objects contained within the system.
        """

        # Since system variables are being set, all analysis objects must be recomputed
        self.flume_sys.reset_analysis_flags(it_counter)

        # Loop through the design variables for the system and set for their components
        for var in self.flume_sys.design_vars_info:
            # Get the indices for the current variable
            start = self.indices[var]["start"]
            end = self.indices[var]["end"]

            # Extract the variable value
            if end - start == 1:
                x_i = x[start:end].item()
            else:
                x_i = x[start:end]

            # ic(x_i)

            # Extract the local name of the variable
            local_name = self.flume_sys.design_vars_info[var]["local_name"]

            # Declare the variable as a design variable for Flume
            self.flume_sys.design_vars_info[var]["instance"].declare_design_vars(
                variables=[local_name]
            )

            # Set the variable value for the Flume analysis object
            self.flume_sys.design_vars_info[var]["instance"].set_var_values(
                variables={local_name: x_i}
            )

        return

    def evalObjCon(self, x):
        """
        Evaluates the objective and constraints for the problem using the current design variables.
        """

        # print("\nCALLING EVALOBJCON")

        # Call the update function, if it was provided
        if self.update is not None:
            self.update(it_num=self.it_counter)

        # If the Flume system does not have an FOI attribute, set it
        if not hasattr(self.flume_sys, "foi"):
            self.flume_sys.declare_foi(global_foi_name=[])

        # TODO: need to add a copy of the top-level analysis list here and remove the analyses after they are performed

        # Check to make sure that the objective analysis info has been set, otherwise raise an error
        if not hasattr(self.flume_sys, "obj_analysis"):
            raise RuntimeError(
                f"The objective information for the system named '{self.flume_sys.sys_name}' has not yet been declared, so evalObjCon can not be executed. Ensure that the function 'declare_objective' has been called."
            )

        # Set the variable values for the various analyses
        self.set_system_variables(x, self.it_counter)

        # Perform the analysis for the objective function
        self.flume_sys.obj_analysis.analyze()

        # Extract the objective function output
        self.obj_name = self.flume_sys.obj_local_name
        # ic(self.flume_sys.obj_analysis.outputs)
        obj = (
            self.flume_sys.obj_analysis.outputs[self.obj_name].value
            * self.flume_sys.obj_scale
        )

        # Loop through each of the constraints and perform their respective analyses
        con_list = []
        for con in self.flume_sys.con_info:
            # ic(con)
            # Perform the analysis for the current constraint function
            self.flume_sys.con_info[con]["instance"].analyze()

            # Extract the output for the constraint
            con_name = self.flume_sys.con_info[con]["local_name"]

            con_val = self.flume_sys.con_info[con]["instance"].outputs[con_name].value

            # ic(con_val)

            # Using the direction information about the constraint, set the normalized constraint value
            if self.flume_sys.con_info[con]["direction"] == "geq":
                # If rhs is not 0.0, scale the constraint
                if self.flume_sys.con_info[con]["rhs"] != 0.0:
                    rhs_val = self.flume_sys.con_info[con]["rhs"]
                    con_val = con_val / rhs_val - 1.0

            elif self.flume_sys.con_info[con]["direction"] == "leq":
                # If rhs is not 0.0, scale the constraint
                if self.flume_sys.con_info[con]["rhs"] != 0.0:
                    rhs_val = self.flume_sys.con_info[con]["rhs"]
                    con_val = 1.0 - con_val / rhs_val
                else:
                    # This step is necessary only for 'leq' to convert constraint to proper form, c(x) >= 0.0
                    con_val *= -1.0

            else:
                raise RuntimeError("Constraint direction must be 'geq' or 'leq'.")

            # con_val -= self.flume_sys.con_info[con]["rhs"]

            # Append the constraint value to the constraints list
            con_list.append(con_val)

        # Call the logger function
        self.flume_sys.log_information(iter_number=self.it_counter)

        # Call the callback function, if it was provided
        if self.callback is not None:
            self.callback(x, self.it_counter)

        # Set the failure flag
        fail = 0

        # Update the iteration counter
        self.it_counter += 1

        return fail, obj, con_list

    def evalObjConGradient(self, x, g, A):
        """
        Evaluates the objective and constraint gradients for the system.
        """

        # print("\nCALLING EVALOBJCONGRADIENT")

        # Check to make sure that the objective analysis info has been set, otherwise raise an error
        if not hasattr(self.flume_sys, "obj_analysis"):
            raise RuntimeError(
                f"The objective information for the system named '{self.flume_sys.sys_name}' has not yet been declared, so evalObjCon can not be executed. Ensure that the function 'declare_objective' has been called."
            )

        # # Set the variable values for the various analyses
        # self.set_system_variables(x)

        # Compute the gradient of the objective function, where the seed value is set to 1.0 for the output of interest
        self.flume_sys.obj_analysis._add_output_seed(outputs=[self.obj_name], seed=1.0)

        self.flume_sys.obj_analysis.analyze_adjoint(debug_print=False)

        # Extract the derivative of the objective wrt the design variables
        for var in self.flume_sys.design_vars_info:
            # Get the indices for the current variable
            start = self.indices[var]["start"]
            end = self.indices[var]["end"]

            # Extract the local name of the variable
            local_name = self.flume_sys.design_vars_info[var]["local_name"]

            # Extract the derivative for the current design variable
            gradx_i = (
                self.flume_sys.design_vars_info[var]["instance"]
                .variables[local_name]
                .deriv
            )

            # Assign the gradient
            g[start:end] = gradx_i * self.flume_sys.obj_scale

        con_index = 0
        # Loop through the constraints in the system
        for con in self.flume_sys.con_info:
            # Extract the local name of the constraint
            con_name = self.flume_sys.con_info[con]["local_name"]
            # ic(con_name)

            # Set the seed for the current constraint
            con_val = self.flume_sys.con_info[con]["instance"].outputs[con_name].value

            if isinstance(con_val, float):
                seed = 1.0
            elif isinstance(con_val, np.ndarray):
                seed = np.ones_like(con_val)
            else:
                raise RuntimeError(
                    f"The type for the constraint value '{con_name}' is not a float or NumPy array, which is unexpected behavior."
                )

            # Add the output seed
            self.flume_sys.con_info[con]["instance"]._add_output_seed(
                outputs=[con_name], seed=seed
            )

            # Perform the adjoint analysis
            self.flume_sys.con_info[con]["instance"].analyze_adjoint(debug_print=False)

            # Loop through the variables in the system
            for var in self.flume_sys.design_vars_info:
                # Get the indices for the current variable
                start = self.indices[var]["start"]
                end = self.indices[var]["end"]

                # Extract the derivative value for the current constraint and variable combination
                local_var_name = self.flume_sys.design_vars_info[var]["local_name"]
                gradc_i = (
                    self.flume_sys.design_vars_info[var]["instance"]
                    .variables[local_var_name]
                    .deriv
                )

                # ic(gradc_i)

                # Assign the constraint gradient, accounting for the scaling, when necessary
                if self.flume_sys.con_info[con]["direction"] == "geq":
                    # If rhs is not 0.0, apply the scale to the constraint
                    if self.flume_sys.con_info[con]["rhs"] != 0.0:
                        rhs_val = self.flume_sys.con_info[con]["rhs"]
                        gradc_i /= rhs_val

                elif self.flume_sys.con_info[con]["direction"] == "leq":
                    # If rhs is not 0.0, scale the constraint
                    if self.flume_sys.con_info[con]["rhs"] != 0.0:
                        rhs_val = self.flume_sys.con_info[con]["rhs"]
                        gradc_i /= -rhs_val
                    else:
                        gradc_i *= -1.0

                A[con_index][start:end] = gradc_i

            # Update the constraint index
            con_index += 1

        # Add the profiling information for the current iteration
        self.flume_sys.profile_iteration(self.it_counter - 1)

        return 0

    def construct_paropt_problem(self) -> ParOptProb:
        """
        Function that creates the ParOpt problem that will be used for optimization
        """

        # Construct the ParOptProblem for the Flume system
        paroptprob = ParOptProb(MPI.COMM_SELF, prob=self)

        return paroptprob

    def get_paropt_default_options(self, output_prefix, algorithm="tr", maxit=1000):
        """
        Get the default options for paropt.
        """

        options = {
            "algorithm": algorithm,
            "tr_init_size": 0.05,
            "tr_min_size": 1e-6,
            "tr_max_size": 10.0,
            "tr_eta": 0.25,
            "tr_infeas_tol": 1e-6,
            "tr_l1_tol": 1e-3,
            "tr_linfty_tol": 0.0,
            "tr_adaptive_gamma_update": True,
            "tr_max_iterations": maxit,
            "mma_max_iterations": maxit,
            "mma_init_asymptote_offset": 0.2,
            "max_major_iters": 100,
            "penalty_gamma": 1e3,
            "qn_subspace_size": 10,
            "qn_type": "bfgs",
            "abs_res_tol": 1e-8,
            "starting_point_strategy": "affine_step",
            # "barrier_strategy": "mehrotra_predictor_corrector",
            "barrier_strategy": "mehrotra",
            "use_line_search": False,
            # "mma_constraints_delta": True,
            "mma_move_limit": 0.1,
            "output_file": os.path.join(output_prefix, "paropt.out"),
            "tr_output_file": os.path.join(output_prefix, "paropt.tr"),
            "mma_output_file": os.path.join(output_prefix, "paropt.mma"),
        }

        return options
