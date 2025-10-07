from scipy.optimize import minimize, Bounds, NonlinearConstraint, BFGS, SR1
import numpy as np
from flume.base_classes.system import System
from icecream import ic


class FlumeScipyInterface:

    def __init__(self, flume_sys: System, callback=None, update=None):
        """
        DOCS:
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
            # Get the keys for the constraints
            con_keys = self.flume_sys.con_info.keys()

            ncon = 0
            for con in con_keys:
                if self.flume_sys.con_info[con]["direction"] == "both":
                    ncon += 2
                else:
                    ncon += 1

            self.ncon = ncon
        else:
            self.ncon = 0

        # Set the number of design variables
        self._set_ndvs()

        return

    def _set_ndvs(self):
        """
        This needs to put the design variables into a flattened list, as SciPy expects that the design variables are flattened into a 1D array that is passed to value/gradient functions.
        """

        # Initialize a dictionary that will store the variable sizes
        self.indices = {}
        ndvs = 0

        # Loop through each of the variables in the System vars dictionary, assembling/flattening into an array
        for var in self.flume_sys.design_vars_info:
            # Extract the instance of the object for the current design variable
            instance = self.flume_sys.design_vars_info[var]["instance"]

            # Get the local name for the variable
            local_name = self.flume_sys.design_vars_info[var]["local_name"]

            # Get the value of the variable
            var_value = instance.variables[local_name].value

            # Extract the size info and store in the variable sizes dictionary
            if isinstance(var_value, float):
                ndvs_i = 1

            elif isinstance(var_value, np.ndarray):
                # Set the size to the size of the NumPy array
                ndvs_i = var_value.size

            else:
                raise RuntimeError(
                    "Variable value is neither a float or NumPy array, which is unexpected behavior."
                )

            # Store the start and end indices for the current design variable, which is used to set/extract the variables from the global design variable vector for SciPy
            self.indices[var] = {"start": ndvs, "end": ndvs + ndvs_i}

            # Add the number of dvs to the total number od design variables
            ndvs += ndvs_i

        # Set the attribute for the number of dvs
        self.ndvs = ndvs

        return

    def _get_dv_bounds(self, x):
        """
        Sets the bounds for the design variables contained in x based on the information associated with the Flume system. Returns a SciPy optimize 'Bounds' object.
        """

        # Initialize the lb and ub objects, which by default have no lower/upper bound
        lb = np.ones_like(x) * -np.inf
        ub = np.ones_like(x) * np.inf

        # Loop through the design variables in the system and set the bounds
        for var in self.flume_sys.design_vars_info:
            # Extract the local name for the current variable
            local_name = self.flume_sys.design_vars_info[var]["local_name"]

            # Extract the start/end indices for the current variable in the global design variable array
            start = self.indices[var]["start"]
            end = self.indices[var]["end"]

            # Set the lower bound value for the current variable, if it exists
            if "lb" in self.flume_sys.design_vars_info[var]:
                lb[start:end] = self.flume_sys.design_vars_info[var]["lb"]
            else:
                pass

            # Set the upper bound value for the current variable, if it exists
            if "ub" in self.flume_sys.design_vars_info[var]:
                ub[start:end] = self.flume_sys.design_vars_info[var]["ub"]
            else:
                pass

        # Construct the Bounds object
        bounds = Bounds(lb=lb, ub=ub)

        return bounds

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

    def _objective_func(self, x, method):
        """
        Computes the objective function value for the system using the current values for the design variables, x.
        """

        # Call the update function, if it was provided
        if self.update is not None:
            self.update(it_num=self.it_counter)

        # Check to make sure that the objective analysis info has been set, otherwise raise an error
        if not hasattr(self.flume_sys, "obj_analysis"):
            raise RuntimeError(
                f"The objective information for the system named '{self.flume_sys.sys_name}' has not yet been declared, so evalObjCon can not be executed. Ensure that the function 'declare_objective' has been called."
            )

        # Set the variable values for the various analyses
        self.set_system_variables(x, self.it_counter)

        # Perform the analysis for the objective function
        self.flume_sys.obj_analysis.analyze(debug_print=False)

        # Extract the objective function output
        self.obj_name = self.flume_sys.obj_local_name
        obj = (
            self.flume_sys.obj_analysis.outputs[self.obj_name].value
            * self.flume_sys.obj_scale
        )

        # Update the iteration counter
        self.it_counter += 1

        if method == "trust-constr":
            print("hello")
            _ = self._constraint_funcs(x, method)

        return obj

    def _grad_objective_func(self, x, method):
        """
        Evaluates the gradient of the objective function for the system
        """

        # Check to make sure that the objective analysis info has been set, otherwise raise an error
        if not hasattr(self.flume_sys, "obj_analysis"):
            raise RuntimeError(
                f"The objective information for the system named '{self.flume_sys.sys_name}' has not yet been declared, so evalObjCon can not be executed. Ensure that the function 'declare_objective' has been called."
            )

        # Compute the gradient of the objective function, where the seed value is set to 1.0 for the output of interest
        self.flume_sys.obj_analysis._add_output_seed(outputs=[self.obj_name], seed=1.0)

        self.flume_sys.obj_analysis.analyze_adjoint(debug_print=False)

        # Initialize the gradient vector according to the number of design vars
        g = np.zeros(self.ndvs)

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

        return g

    def _constraint_funcs(self, x, method):
        """
        Computes the constraint function values for the system using the current values for the design variables, x.
        """

        # Detect the method used for the optimization
        if method in ["SLSQP"]:
            # Set the boolean flag for detecting the method type
            dict_method = True

        elif method in ["COBYQA", "trust-constr", "COBYLA"]:
            # Set the boolean flag for detecting the method type
            dict_method = False

            if method == "trust-constr":
                raise NotImplementedError(
                    f"Method '{method}' has not been implemented yet."
                )
        else:
            raise RuntimeError(
                f"The provided method '{method}' is not supported with constraints with SciPy optimize minimize. Valid options are 'COBYLA', 'COBYQA', 'SLSQP', and 'trust-constr' (note trust-constr not implemented yet)."
            )

        # Loop through each of the constraints and perform their respective analyses
        con_list = []
        for con in self.flume_sys.con_info:
            # Extract the sub-dictionary of information associated with the current constraint
            con_i_info = self.flume_sys.con_info[con]

            # Get the constraint type for SciPy
            if con_i_info["direction"] == "geq" or con_i_info["direction"] == "leq":
                con_type = "ineq"
            elif con_i_info["direciton"] == "both":
                con_type = "eq"
            else:
                raise RuntimeError(
                    f"The value for 'direction' for the constraint {con} must be 'geq', 'leq', or 'both' and not {con_i_info['direction']}"
                )

            # Construct a lambda function for the constraint wrapper
            con_func = lambda x: self._constraint_wrapper(x, con_i_info)

            # Construct a lambda function for the constraint Jacobian wrapper
            con_jac_func = lambda x: self._constraint_jac_wrapper(x, con_i_info)

            # Add the constraint
            if dict_method:
                # Construct the dictionary for the current constraint
                con_dict = {"type": con_type, "fun": con_func, "jac": con_jac_func}

                # Append to the list
                con_list.append(con_dict)
            else:
                # Construct the NonlinearConstraint object
                con_obj = NonlinearConstraint(
                    fun=con_func, lb=0.0, ub=np.inf, jac=con_jac_func
                )

                # Append to the list
                con_list.append(con_obj)

        # Return the constraint list
        return con_list

    def _constraint_wrapper(self, x, con_i_info: dict):
        """
        DOCS:
        """

        # Extract the instance, direction, and rhs value associated with the constraint
        instance = con_i_info["instance"]
        direction = con_i_info["direction"]
        rhs = con_i_info["rhs"]
        local_name = con_i_info["local_name"]

        # Perform the analysis for the current constraint
        instance.analyze(debug_print=False)

        # Extract the constraint value
        con_val = instance.outputs[local_name].value

        # Using the direction information about the constraint, set the normalized constraint value
        if direction == "geq":
            # If rhs is not 0.0, scale the constraint
            if rhs != 0.0:
                con_val = con_val / rhs - 1.0

            # If the rhs is < 0, flip the sign for the constraint
            if rhs < 0.0:
                con_val *= -1.0

        elif direction == "leq":
            # If rhs is not 0.0, scale the constraint
            if rhs != 0.0:
                con_val = 1.0 - con_val / rhs
            else:
                # This step is necessary only for 'leq' to convert constraint to proper form, c(x) >= 0.0
                con_val *= -1.0

            # If the rhs is < 0, flip the sign for the constraint
            if rhs < 0.0:
                con_val *= -1.0

        elif direction == "both":
            # If rhs is not 0.0, scale the constraints
            if rhs != 0.0:
                # Normalized constraint value
                con_val = con_val / rhs - 1.0

                # If the rhs is < 0, flip the sign for the constraint
                if rhs < 0.0:
                    con_val *= -1.0

            else:
                # Unormalized constraint value
                con_val = con_val

        else:
            raise RuntimeError("Constraint direction must be 'geq', 'leq', or 'both'.")

        return con_val

    def _constraint_jac_wrapper(self, x, con_i_info: dict):
        """
        DOCS:
        """

        # Extract the instance, direction, and rhs value associated with the constraint
        instance = con_i_info["instance"]
        direction = con_i_info["direction"]
        rhs = con_i_info["rhs"]
        local_name = con_i_info["local_name"]

        # Get the value for the constraint
        con_val = instance.outputs[local_name].value

        # Set the seed for the current constraint
        if isinstance(con_val, float):
            seed = 1.0
            con_size = 1
        elif isinstance(con_val, np.ndarray):
            seed = np.ones_like(con_val)
            con_size = con_val.size
        else:
            raise RuntimeError(
                f"The type for the constraint value '{local_name}' is not a float or NumPy array, which is unexpected behavior."
            )

        # Add the output seed
        instance._add_output_seed(outputs=[local_name], seed=seed)

        # Perform the adjoint analysis
        instance.analyze_adjoint(debug_print=False)

        # Initialize the constraint Jacobian
        jac = np.zeros((con_size, x.size))

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

            # Modify the constraint gradient, accounting for scaling, when necessary
            if direction == "geq":
                # If rhs is not 0.0, scale the constraint
                if rhs != 0.0:
                    gradc_i /= rhs

                # If the rhs is < 0, flip the sign
                if rhs < 0.0:
                    gradc_i *= -1.0

            elif direction == "leq":
                # If rhs is not 0.0, scale the constraint
                if rhs != 0.0:
                    gradc_i /= -rhs
                else:
                    # This step is necessary only for 'leq' to convert the constraint to the proper form, c(x) >= 0.0
                    gradc_i *= -1.0

                # If the rhs is < 0, flip the sign for the constraint
                if rhs < 0.0:
                    gradc_i *= -1.0

            elif direction == "both":
                # If rhs is not 0.0, scale the constraint
                if rhs != 0.0:
                    gradc_i /= rhs

                # If the rhs is < 0, flip the sign for the constraint
                if rhs < 0.0:
                    gradc_i *= -1.0

            else:
                raise RuntimeError(
                    "Constraint direction must be 'geq', 'leq', or 'both'."
                )

            # Assign the value of gradc_i to the Jacobian
            jac[:, start:end] = gradc_i

        return jac

    def get_default_options(self, maxit=100):
        """
        Get the dictionary of default options that are used when calling SciPy minimize.
        """

        # Construct the dictionary with the default options
        options = {"disp": True, "maxiter": maxit, "ftol": 1e-10, "verbose": 2}

        return options

    def optimize_system(self, x0, options=None, method="SLSQP", maxit=100):
        """
        Function that the user calls to execute SciPy's optimize minimize function.
        """

        # Get the default options if the user did not provide any
        if options is not None:
            options = self.get_default_options(maxit=maxit)
            pass

        # Call the minimize function for the system
        res = minimize(
            fun=self._objective_func,
            x0=x0,
            jac=self._grad_objective_func,
            args=method,
            bounds=self._get_dv_bounds(x0),
            constraints=self._constraint_funcs(x0, method),
            options=options,
            method=method,
        )

        return res.x, res
