# Component sourced from OpenMDAO docs on 09/08/2026

import openmdao.api as om


class Paraboloid(om.ExplicitComponent):
    """
    Evaluates the equation f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3.
    """

    def setup(self):
        self.add_input("x", val=0.0, desc="x value")
        self.add_input("y", val=0.0, desc="y value")

        self.add_output("f_xy", val=0.0)

    def setup_partials(self):
        # Finite difference all partials.
        self.declare_partials("*", "*", method="fd")

    def compute(self, inputs, outputs):
        """
        f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3

        Minimum at: x = 6.6667; y = -7.3333
        """
        x = inputs["x"]
        y = inputs["y"]

        outputs["f_xy"] = (x - 3.0) ** 2 + x * y + (y + 4.0) ** 2 - 3.0
