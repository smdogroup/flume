import numpy as np
import openmdao.api as om
from flume.base_classes.analysis import Analysis
from flume.base_classes.state import State
from icecream import ic


class OpenMDAOComponentAnalysis(Analysis):

    def __init__(
        self,
        om_component: om.ExplicitComponent,
        obj_name: str,
        sub_analyses: list = [],
        **kwargs,
    ):

        # Default parameters
        self.default_parameters = {}

        # Store the OpenMDAO component
        self.om_component = om_component

        # Construct an OpenMDAO problem that solely contains the input component in the model
        self.prob = om.Problem()
        self.prob.model.add_subsystem(obj_name, subsys=self.om_component)

        # Perform the setup for the problem so that the component's inputs and outputs are accessibl
        self.prob.setup()
        self.prob.final_setup()

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Get the input and output names
        om_inputs = self.prob.model.list_inputs(out_stream=None, val=True, desc=True)
        om_outputs = self.prob.model.list_outputs(out_stream=None, val=True, desc=True)

        # Initialize the dictionary that maps Flume names to OpenMDAO names
        self._var_map = {}
        self._var_rev_map = {}

        # Initialize the list that stores all OpenMDAO variable names
        self._om_var_names = []

        # Setup the variables dictionary for the Flume object (maps local Flume names to OpenMDAO model name)
        self.variables = {}

        for var in om_inputs:
            # Get the variable name and variable info from the OpenMDAO Component
            var_name = var[0]
            var_info = var[1]

            # Append the variable name to the list
            self._om_var_names.append(var_name)

            # Get the local variable name (split from global_name.local_name)
            local_name = var_name.split(".")[-1]

            value = self._to_native(var_info["val"])
            desc = var_info["desc"]

            # Add the entry to the var map
            self._var_map[local_name] = var_name

            # Construct the State object
            var_state = State(value=value, desc=desc, source=self)

            # Add the State object to the variables dictionary
            self.variables[local_name] = var_state

            # Add the entry to the reverse map (OpenMDAO -> Flume)
            self._var_rev_map[var_name] = local_name

        # Setup the dictionary that maps Flume output names to OpenMDAO output names
        self._out_map = {}
        self._out_rev_map = {}

        # Initialize the list that stores all OpenMDAO output names
        self._om_out_names = []

        # Loop through the outputs and populate the output map dictionary and list
        for out in om_outputs:
            # Get the OpenMDAO output name
            out_name = out[0]
            out_info = out[1]

            # Append the output name to the list
            self._om_out_names.append(out_name)

            # Get the local output name
            local_name = out_name.split(".")[-1]

            # Add the entry to the output name map
            self._out_map[local_name] = {}
            self._out_map[local_name]["om_name"] = out_name
            self._out_map[local_name]["desc"] = out_info["desc"]

            # Add the entry to the reverse map (OpenMDAO -> Flume)
            self._out_rev_map[out_name] = local_name

        return

    @staticmethod
    def _to_native(val):
        """
        Convert an OpenMDAO value into Flume's native representation.

        OpenMDAO stores every variable/output as a NumPy array (a scalar is a shape-(1,) array). Flume's convention, however, is that scalar quantities are Python scalars (e.g. States are built with value=0.0), and the optimizer interfaces rely on this..

        This helper returns a Python scalar for size-1 quantities and an independent NumPy array copy otherwise. Returning a copy (or an immutable scalar) is also required because prob.get_val hands back a live reference into OpenMDAO's internal vectors, which are overwritten in place on the next run_model.
        """
        arr = np.asarray(val)
        if arr.size == 1:
            # .item() preserves the scalar's own type (float or complex)
            return arr.item()
        return np.array(arr)

    def _analyze(self):
        """
        Private analyze method, which calls run_model method to compute the Component outputs using the variable (input) values, and then assigns them into the States.
        """

        # Loop through the variables and set the values in the State object in the OpenMDAO model
        for var in self.variables:
            # Get the value of the variable State
            var_val = self.variables[var].value

            # Get the name of the OpenMDAO variable
            om_name = self._var_map[var]

            # Set the variable value in the OpenMDAO model
            self.prob.set_val(name=om_name, val=var_val)

        # Execute run_model for the problem
        self.prob.run_model()

        # Extract the outputs and assign in the outputs dictionary
        self.outputs = {}

        for out in self._out_map:
            # Get the OpenMDAO output name
            om_out_name = self._out_map[out]["om_name"]
            desc = self._out_map[out]["desc"]

            # Get the value for the OpenMDAO output, converting to Flume's native
            # representation
            out_val = self._to_native(self.prob.get_val(name=om_out_name))

            # Construct the State object and assign it in the dictionary
            self.outputs[out] = State(value=out_val, desc=desc, source=self)

        return

    def _analyze_adjoint(self):
        """
        Private adjoint Analysis, which calls the compute_totals method to compute the total derivatives of each output wrt each input, and then accumulates across outputs for updating the variables' derivatives.
        """

        # Execute the compute_totals method for the OpenMDAO model to compute the total derivatives through the Component
        derivs = self.prob.compute_totals(of=self._om_out_names, wrt=self._om_var_names)

        # Loop through all outputs
        for out in self.outputs:
            # Get the adjoint value for the output
            outb = self.outputs[out].deriv

            # Get the OpenMDAO output name
            om_out_name = self._out_map[out]["om_name"]

            # Loop through all variables and accumulate the contributions to each variable's derivative
            for var in self.variables:

                # Get the existing derivative value
                varb = self.variables[var].deriv

                # Get the partial derivative for the current (output, variable) combination
                om_var_name = self._var_map[var]

                dout_dvar = derivs[
                    (om_out_name, om_var_name)
                ]  # Note this is a 2d Jacobian shaped (n_out, n_var)

                # Reverse-mode accumulation x_bar += J^T @ f_bar
                J = np.atleast_2d(dout_dvar)  # (n_out, n_var)
                fbar = np.atleast_1d(outb)  # (n_out,)
                contrib = J.T @ fbar  # (n_var,)

                # Match the variable's native representation (scalar vs array)
                if not isinstance(self.variables[var].value, np.ndarray):
                    contrib = contrib.item()

                varb = varb + contrib

                # Assign into the dictionary
                self.variables[var].set_deriv_value(varb)

        return
