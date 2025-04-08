from flume.base_classes.analysis import Analysis
import graphviz as gv
from icecream import ic
from typing import List
from flume.interfaces.utils import Logger
import os


class System:
    """
    This is a potential new class that is being workshopped. The idea here is to have this combine one or several analysis objects into a new class. This will help create a wrapper around separate analysis procedures, which could then be used within an optimization problem. Some of the methods associated with this class might include:

    1. Network visualization: construct a plot for visualizing the network (using networkX Python package)
    2. Adding/declaring design variables
    3. Declaring objective/constraint functions for different outputs
    4. Declaring other FOI
    5. Adding FOI, obj, cons to a log file (flume.log) that can be updated throughout the optimization procedure.
    """

    def __init__(
        self,
        sys_name: str,
        top_level_analysis_list: List[Analysis],
        log_name: str = "flume.log",
        log_prefix: str = None,
    ):

        # Store the name for the system
        self.sys_name = sys_name

        # Store the list of analyses
        self.top_level_analysis_list = top_level_analysis_list

        # Store the information for the logging
        self.log_name = log_name
        self.log_prefix = log_prefix

        # Configure the file path for th elog file
        Logger.set_log_path(os.path.join(self.log_prefix, self.log_name))

        # Assemble the list of all analysis objects
        self.full_analysis_list = self.assemble_full_analysis_list()

        # Initialize the constraints dictionary (default assumes no constraints)
        self.con_info = {}

        return

    def reset_analysis_flags(self):
        """
        Resets all of the analysis flags for all analyses within the system to be False. This is done when variable values are updated, as all systems must be analyzed again to propagate chanes in design variable values.
        """

        # Loop through each top-level analysis, resetting each analysis in the stack
        for top_level in self.top_level_analysis_list:

            # If the stack does not exist, make it
            if not hasattr(top_level, "stack"):
                # Assemble the stack
                top_level.stack = top_level._make_stack()

            # Loop through each object within the current top-level's analysis stack and set the analyzed attribute to be False
            for analysis in top_level.stack:
                analysis.analyzed = False

        return

    def assemble_full_analysis_list(self):
        """
        Assembles the full list of analyses that comprise the overall system architecture.
        """

        # Initialize the analysis list
        full_analysis_list = []

        # Loop through each top-level analysis, assemble the stack, and append the analyses to the total analysis list
        for top_level in self.top_level_analysis_list:
            # Assemble the stack
            stack = top_level._make_stack()

            # Add the entries in the current stack to the list if they are not there already
            for analysis in stack:
                if analysis not in full_analysis_list:
                    full_analysis_list.append(analysis)
                else:
                    continue

        return full_analysis_list

    def graph_network(self, make_connections=False):
        """
        DOCS: goal here is to construct the directed acyclic graph using the information for each of the entries in the analyses_list

        make_connections will likely have to call analyze on the system to construct the connections object (since outputs do not exist by default)
        """

        # Create the graph
        # graph = nx.Graph()
        graph = gv.Digraph(
            name=f"{self.sys_name.upper()}", graph_attr={"rankdir": "LR"}
        )

        # Initialize an empty list to store the systems added to the graph as nodes
        self.nodes = []

        # Loop through and add nodes to the graph (outer loop is top-level analyses, inner loop is individual sub-analyses)
        for analysis in self.top_level_analysis_list:
            # Check if the object has a stack already, otherwise assemble
            if hasattr(analysis, "stack"):
                stack = analysis.stack
                pass
            else:
                stack = analysis._make_stack()

            # Check if the object has been connected already, otherwise perform the analysis to establish connections map
            if analysis.connected:
                pass
            else:
                analysis.analyze()

            # Loop through sub analyses in the current stack and add the nodes if they do not already exist
            for i, sub in enumerate(stack):
                # Add node if it is not in the system network already
                if sub not in self.nodes:
                    self.nodes.append(sub)

                    # Set the color depending on whether the analysis is top-level
                    if sub in self.top_level_analysis_list:
                        # Extract the names of the outputs
                        out_labels = list(sub.outputs.keys())
                        out_str = ", ".join(out_labels)

                        # Add the node
                        graph.node(
                            sub.obj_name,
                            f"{sub.obj_name}\nOutputs: {out_str}",
                            color="red",
                        )
                    else:
                        graph.node(sub.obj_name, f"{sub.obj_name}")
                else:
                    pass

                # Add edge if not the first entry in the stack
                if i == 0:
                    pass
                else:
                    # Extract the States that are connected between the objects
                    connect_labels = list(sub.connects.keys())
                    connect_str = ", ".join(connect_labels)

                    # Graph the edge and add the label for the connection
                    graph.edge(
                        stack[i - 1].obj_name,
                        sub.obj_name,
                        label=f"{connect_str}",
                    )

        return graph

    def _find_analysis_object(self, instance_name, var_name):
        """
        Searches through the full analysis list for the system and finds the instance of the analysis object associated with the argument for the instance_name.
        """

        # Loop through the full analysis list
        for analysis in self.full_analysis_list:
            # Check if the current object name matches the provided instance name and return if so
            if analysis.obj_name == instance_name:
                return analysis

        # Raise a RuntimeError if this point is reached, as the instance name did not find a match
        raise RuntimeError(
            f"No instance found for object named '{instance_name}'! Verify the definition for the State named '{instance_name}.{var_name}'"
        )

    def declare_objective(self, global_obj_name):
        """
        Sets the objective function for the optimization problem according to the provided global output name. Nominally, this output should be associated with one of the top-level analyses for the system (i.e. included in the analysis_list).
        """

        # Using the provided objective name, store the associated analysis object and the local variable name
        obj_analysis_name, self.obj_local_name = global_obj_name.split(".")

        # Find the instance for the objective analysis object
        self.obj_analysis = self._find_analysis_object(
            instance_name=obj_analysis_name, var_name=self.obj_local_name
        )

        return

    def declare_constraints(self, global_con_name: dict):
        """
        Sets the constraints for the optimization problem according to the provided global output names. Following the convention for ParOpt, the constraints are assumed to be of the following form:

        c(x) >= 0.0

        Parameters
        ----------
        global_con_name_dict : dict
            This is a dictionary. The keys correspond to the global output names for the constraints. The values correspond to the right hand side for the constraint following the format above. If the value is not 0.0, then the RHS value will be subtracted from the evaluated constraint value during the optimization procedure.

        Example
        -------

        Here, the constraints are defined as follows:
            x >= 0.0
            y >= 1.0

        Thus, the argument here is given as: global_con_name = {"block1.x":0.0, "block2.y":1.0}
        """

        # Loop through the keys in the dictionary and add them to the constraints for the system
        for key in global_con_name.keys():
            # Add the key to the con_info dictionary
            self.con_info[key] = {}

            # Split the string
            con_analysis_name, con_local_name = key.split(".")

            # Find the analysis object for the current constraint
            con_analysis = self._find_analysis_object(
                con_analysis_name, var_name=con_local_name
            )

            # Add the analysis object and local constraint name for the con_info dictionary
            self.con_info[key]["instance"] = con_analysis
            self.con_info[key]["local_name"] = con_local_name
            self.con_info[key]["rhs"] = global_con_name[key]

        return

    def declare_design_vars(self, global_var_name: dict):
        """
        Sets the design variables for the system according to the provided global variable names.

        Parameters
        ----------
        global_var_name : dict
            This is a dictionary of dictionaries. Keys for the first dictionary correspond to the global variable names that should be added to the set of design variables, and the values contain information about the variable bounds. If the inner dictionary is empty, then no bounds are provided. Otherwise, the values for the variables lower bound 'lb' and upper bound 'ub' will be added, if included.

        Example
        -------

        Here, the bounds are defined as follows:
            1.0 <= x <= 2.0
            y has no bounds
            0.0 <= z

        Thus, the argument here is: global_var_name = {"block1.x":{"lb":1.0, "lb":2.0}, "block2.y":{}, "block3.z":{"lb":0.0}}
        """

        # Initialize the design variables dictionary
        self.design_vars_info = {}

        # Loop through the keys in the dictionary and add them to the design variables for the system
        for key in global_var_name:
            # Add the key to the design_vars_info dictionary
            self.design_vars_info[key] = {}

            # Split the string for the current variable
            var_analysis_name, var_local_name = key.split(".")

            # Find the analysis object instance for the current design variable
            var_analysis = self._find_analysis_object(var_analysis_name, var_local_name)

            # Add the analysis object and local variable name to the dictionary
            self.design_vars_info[key]["instance"] = var_analysis
            self.design_vars_info[key]["local_name"] = var_local_name

            # Add the bounds, if they were specified
            if global_var_name[
                key
            ]:  # This is evaluated as True if the variable has at leaast one bound specified

                # Add the lower bound, if specified
                if "lb" in global_var_name[key].keys():
                    self.design_vars_info[key]["lb"] = global_var_name[key]["lb"]

                # Add the upper bound, if specified
                if "ub" in global_var_name[key].keys():
                    self.design_vars_info[key]["ub"] = global_var_name[key]["ub"]
            else:
                continue

        return

    def declare_foi(self, global_foi_name: list):
        """
        Sets the functions of interest that are to be tracked by the logger at each iteration. Here, the names must be provided according to their global names. If already declared, the objective and constraints are added by default, but other outputs can be included.
        """

        # Initialize the foi dictionary
        foi = {}

        # Add the objective if it has already been declared
        if hasattr(self, "obj_analysis"):
            foi["obj"] = {}

            foi["obj"]["instance"] = self.obj_analysis
            foi["obj"]["local_name"] = self.obj_local_name

        else:
            raise RuntimeError("No objective has been declared!")

        # Add the constraints if it has already been declared
        if hasattr(self, "con_info"):
            foi["cons"] = {}

            # For each constraint in self.con_info, add it to the foi dictionary
            for key in self.con_info:
                foi["cons"][key] = {}

                foi["cons"][key]["instance"] = self.con_info[key]["instance"]
                foi["cons"][key]["local_name"] = self.con_info[key]["local_name"]

        # Add the names in the global_foi_name input list to the foi dictionary
        foi["other"] = {}

        for name in global_foi_name:
            # Split the name for the global foi
            foi_analysis_name, foi_local_name = name.split(".")

            # Find the analysis object associated with the global_foi_name
            foi_analysis = self._find_analysis_object(
                instance_name=foi_analysis_name, var_name=foi_local_name
            )

            # Add the analysis object and local state name to the foi dictionary
            foi["other"][name] = {}
            foi["other"][name]["instance"] = foi_analysis
            foi["other"][name]["local_name"] = foi_local_name

        # Store the foi dictionary
        self.foi = foi

        return

    def log_information(self, iter_number):
        """
        DOCS:
        """

        # Log the header names if the current iter number is divisible by 10
        if iter_number % 10 == 0:
            # Log the header for the iter number and objective
            obj_header = f"obj: {self.obj_local_name}"
            Logger.log("\n%5s%20s" % ("iter", obj_header), end="")

            # Log the constraints
            for con in self.foi["cons"].keys():
                con_header = f"con: {self.foi['cons'][con]['local_name']}"
                Logger.log("%20s" % con_header, end="")

            # Log the other functions of interest
            for other in self.foi["other"].keys():
                other_header = f"other: {self.foi['other'][other]['local_name']}"
                Logger.log("%20s" % other_header, end="")

        # Log the values for the current iteration and objective value
        obj_val = (
            self.foi["obj"]["instance"].outputs[self.foi["obj"]["local_name"]].value
        )

        Logger.log("\n%5d%20.10e" % (iter_number, obj_val), end="")

        # Log the values for the constraints at the current iter
        for con in self.foi["cons"].keys():
            con_val = (
                self.foi["cons"][con]["instance"]
                .outputs[self.foi["cons"][con]["local_name"]]
                .value
            )

            if not isinstance(con_val, str):
                con_val = "%20.10e" % con_val

            Logger.log("%20s" % con_val, end="")

        # Log the values for the other FOI
        for other in self.foi["other"].keys():
            other_val = (
                self.foi["other"][other]["instance"]
                .outputs[self.foi["other"][other]["local_name"]]
                .value
            )

            if not isinstance(other_val, str):
                other_val = "%20.10e" % other_val

            Logger.log("%20s" % other_val, end="")

        return
