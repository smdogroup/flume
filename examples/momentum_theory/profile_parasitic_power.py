import numpy as np
from base_classes.analysis_base import Analysis, State


class ProfileParasiticPower(Analysis):
    """
    Analysis class that computes the profile and parasitic power terms noramlized by advance ratio (influence of flight speed handled by momentum theory segment).
    """

    def __init__(self, obj_name: str, **kwargs):
        """
        Analysis class for computing profile and parasitic power terms normalized by their advance ratio contributions.

        Parameters
        ----------
        obj_name : str
            Name for the analysis object

        Keyword Arguments
        -----------------
        R : float
            Rotor radius (meters)
        """

        # Set the default parameters for the profile power analysis
        self.default_parameters = {
            "R": 1.0,  # rotor radius (m),
        }

        # Set the parameters for the segment using the default arguments or the kwargs
        self.parameters = {}
        for key in self.default_parameters:
            if key in kwargs:
                if type(kwargs[key]) == type(self.default_parameters[key]):
                    self.parameters[key] = kwargs[key]
                else:
                    raise TypeError(
                        f"Type for parameter '{key}' does not match the expected type. Input type for {key} was {type(kwargs[key])} but expected is {type(self.default_parameters[key])}."
                    )
            else:
                self.parameters[key] = self.default_parameters[key]

        # Set the default values for the variables, i.e. the inputs that can change
        f_var = State(
            value=10.0,
            desc="Equivalent flat-plate area of the helicopter, which is used to compute the parasitic power (units: m**2). This could be determined via drag synthesis or CFD.",
            source=self,
        )

        sigma_var = State(
            value=0.1,
            desc="Solidity of the rotor blades (assumed rectangular, i.e. constant chord)",
            source=self,
        )

        C_d0_var = State(
            value=0.008,
            desc="Profile drag coefficient of the airfoils that comprise the rotor blades (assumed constant for rectangular blades)",
            source=self,
        )

        self.variables = {
            "f": f_var,
            "sigma": sigma_var,
            "C_d0": C_d0_var,
        }

        # Perform the base class initialization procedure
        super().__init__(obj_name=obj_name, sub_analyses=[], **kwargs)

        return

    def _analyze(self, mode="real"):
        """
        Private analysis method for the object that computes all of the outputs with the provided input variables and parameters to the object.
        """

        if mode == "real":
            dtype = float
        elif mode == "complex":
            dtype = complex
        else:
            raise ValueError("Analysis mode must be either 'real' or 'complex'.")

        # Compute the modifier term for the profile power coefficient
        self.m_P0 = self.variables["sigma"].value * self.variables["C_d0"].value / 8

        # Compute the modifier term for the parasitic power coefficient
        A = np.pi * self.parameters["R"] ** 2

        self.n_Pp = 0.5 * self.variables["f"].value / A

        # Update the analyzed attribute
        self.analyzed = True

        # Assign the outputs to the outputs dictionary
        self.outputs = {}

        self.outputs["m_P0"] = State(
            value=self.m_P0,
            desc="Modifier term for the profile power, which can then be multiplied by the term with the advance ratio to compute the profile power coefficient.",
            source=self,
        )

        self.outputs["n_Pp"] = State(
            value=self.n_Pp,
            desc="Modifier term for the parasitic power, which can then be multiplied by the term with the advance ratio to compute the parasitic power coefficient.",
            source=self,
        )

        return

    def _analyze_adjoint(self):
        """
        Private adjoint analysis method for the object that computes the derivative values for the input variables via the adjoint method.
        """

        # Assert that the adjoint has been initialized
        assert (
            self.adjoint_initialized
        ), f"The adjoint for '{self.__class__.__name__}' named '{self.obj_name}' has not been initialized yet, so the adjoint analysis cannot be performed."

        # Extract the derivative values for the outputs
        m_P0b = self.outputs["m_P0"].deriv
        n_Ppb = self.outputs["n_Pp"].deriv

        # print(f"PPP m_P0b = {self.outputs['m_P0'].deriv}")

        # Compute the derivative values for the variables
        A = np.pi * self.parameters["R"] ** 2

        fb = n_Ppb / (2 * A)

        C_d0b = m_P0b * self.variables["sigma"].value / 8

        sigmab = m_P0b * self.variables["C_d0"].value / 8

        # Set the derivative values for the variables
        self.variables["f"].deriv += fb
        self.variables["C_d0"].deriv += C_d0b
        self.variables["sigma"].deriv += sigmab

        # Update the attribute for the adjoint analysis
        self.adjoint_analyzed = True

        return


if __name__ == "__main__":

    # Create the ProfileParasiticPower object
    PPP = ProfileParasiticPower(obj_name="PPP")

    # Analyze system
    PPP.analyze()

    # Set the design variables
    PPP.declare_design_vars(variables=["f"])

    # Test the adjoint for PPP
    PPP._test_combined_adjoint(debug_print=True)

    # # Get the output values
    # outs = PPP.get_output_values()

    # for key in outs:
    #     print(f"{key} = {outs[key]}")
