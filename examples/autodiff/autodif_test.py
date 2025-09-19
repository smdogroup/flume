from flume.base_classes.analysis import Analysis
from flume.base_classes.state import State
from icecream import ic


class Parabola(Analysis):
    def __init__(self, obj_name: str, sub_analyses=[], **kwargs):
        # Set the default parameters
        self.default_parameters = {"a": 2.0}

        # Perform the base class object initialization
        super().__init__(obj_name=obj_name, sub_analyses=sub_analyses, **kwargs)

        # Set the default states for the variables
        x_var = State(value=1.0, desc="X value", source=self)

        self.variables = {"x": x_var}

        return

    # NOTE: perhaps add an additional function here that is specific for jax; alternatively, perhaps there is a way to add a decorator which would modify the behavior of _analyze, but this might be harder to implement; if the class has the method that is jax compatible, then during analyze adjoint jax's grad function would be called internally, allowing the user to avoid writing the _analyze_adjoint function
    def _analyze(self):

        # Extract variable value
        x = self.variables["x"].value

        # Extract parameter values
        a = self.parameters["a"]

        y = a * x**2

        # Store the outputs
        self.analyzed = True

        self.outputs = {}

        self.outputs["y"] = State(value=y, desc="Y value", source=self)

        return

    def _analyze_adjoint(self):
        return
