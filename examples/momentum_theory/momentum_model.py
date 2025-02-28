import numpy as np
from flume.base_classes.analysis import Analysis, State
import matplotlib.pyplot as plt
import json
from examples.momentum_theory.profile_parasitic_power import ProfileParasiticPower
from typing import Optional


class MomentumTheorySegment(Analysis):
    """
    Analysis class that computes numerous outputs that characterize a segment of a mission profile using momentum theory analysis techniques.
    """

    def __init__(
        self,
        seg_name: str,
        mode: str,
        sub_analysis: Optional[ProfileParasiticPower] = None,
        **kwargs,
    ):
        """
        Analysis class for performing an analysis of a momentum theory segment.

        Parameters
        ----------
        seg_name : str
            Name for the momentum theory segment object
        mode : str
            Operational mode for the segment. Options are one of the following: ['hover', 'axial_climb', 'axial_descent', 'climb', 'descent', 'forward']. This controls the equations used in the analysis
        sub_analysis : None or ProfielParasiticPower object
            Optional argument that specifies whether the profile and parasitic power terms are computed from a sub-analysis

        Keyword Arguments
        -----------------
        mode : str
            Operational mode, inherited from the required input for the class
        Omega : float
            Rotational speed of the rotor blades (rad/s)
        R : float
            Rotor radius (m)
        rho : float
            Density of the atmosphere (kg/m**3)
        weight : float
            Weight of the rotorcraft (N)
        kappa : float
            Induced power factor for the momentum theory analysis
        rc_type : str
            Either 'single' or 'coaxial', which specifies whether the rotorcraft type is a traditional (single rotor disk) rotorcraft or a coaxial rotorcraft
        kappa_int : float
            Interference factor used for the coaxial rotor system
        equilibrium : bool
            Boolean value that dictates whether the forces of flight are in equilibrium for the rotorcraft (if True, rotor disk angle of attack is solved for directly)
        """

        # Store the sub_analysis depending on whether the argument is a class or None
        if sub_analysis is None:
            sub_analyses = []
        elif isinstance(sub_analysis, ProfileParasiticPower):
            sub_analyses = [sub_analysis]
        else:
            raise TypeError(
                "Input for 'sub_analysis' must be either 'None' or an instance of the 'ProfileParasiticPower' class."
            )

        # Set a list of valid mode options
        mode_opts = [
            "hover",
            "axial_climb",
            "climb",
            "axial_descent",
            "descent",
            "forward",
        ]

        # Check to make sure the input mode is a valid option
        if mode not in mode_opts:
            raise ValueError(
                f"Input for mode of '{mode}' is not a valid choice. Valid input modes are {[mode_type for mode_type in mode_opts]}"
            )
        # Set the mode as an attribute if so; this is used to alter calculations if necessary throughout the analysis
        else:
            self.mode = mode

        # Set the default parameters for the segment (note that self.mode is a required argument but is set into parameters dictionary)
        self.default_parameters = {
            "mode": self.mode,  # operational mode
            "Omega": 40.0,  # rotational speed of rotor blades (rad/s)
            "R": 1.0,  # rotor radius (m)
            "rho": 1.293,  # density (kg/m^3)
            "weight": 20.0,  # weight of the rotorcraft (N)
            "kappa": 1.15,  # induced power factor
            "rc_type": "single",  # type for the rotorcraft (coaxial or single)
            "kappa_int": 1.281,  # interference factor used for coaxial rotor system
            "equilibrium": False,  # boolean that dictates whether or not the forces of flight are in equilibrium; if True, rotor disk angle of attack is solved for
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

        # Check to make sure that the rc type is either single or coaxial
        if self.parameters["rc_type"] not in ["single", "coaxial"]:
            raise ValueError(
                f"Input parameter for rotorcraft type (rc_type) of '{self.parameters['rc_type']}' is not a valid option. This parameter must be either 'single' or 'coaxial'."
            )

        # Perform the base class initialization procedure
        super().__init__(obj_name=seg_name, sub_analyses=sub_analyses, **kwargs)

        # Set the variables using the State class for the segment analysis
        alpha_var = State(
            value=0.0, desc="Rotor disk angle of attack in radians", source=self
        )

        V_inf_var = State(value=0.0, desc="Flight speed of the rotorcraft in m/s")

        m_P0_var = State(
            value=0.0,
            desc="Modifier term for the profile power, which can then be multiplied by the term with the advance ratio to compute the profile power coefficient.",
            source=self,
        )

        n_Pp_var = State(
            value=0.0,
            desc="Modifier term for the parasitic power, which can then be multiplied by the term with the advance ratio to compute the parasitic power coefficient.",
            source=self,
        )

        # Set the default State objects for the variables
        self.variables = {
            "alpha": alpha_var,
            "V_inf": V_inf_var,
            "m_P0": m_P0_var,
            "n_Pp": n_Pp_var,
        }

        return

    def _initialize_analysis(self, mode="real"):
        """
        Initializes the internal data for the analysis procedure contained within the object. Also, sets the weight attribute, which is how much weight a single rotor disk carries, depending on the rc_type parameter.
        """

        # If the rotorcraft type is coaxial, set the weight attribute to be half the weight of the rotorcraft (assuming thrust-balanced rotors for the system); otherwise, the weight is unchanged since the single rotor disk must carry all the weight
        if self.parameters["rc_type"] == "coaxial":
            self.weight = self.parameters["weight"] / 2.0
        else:
            self.weight = self.parameters["weight"]

        # Initialize the flag that specifies whether the analysis procedure has been performed
        self.analyzed = False

        return

    def _analyze(self, mode="real"):
        """
        Private analysis method for the object that computes all of the outputs with the provided input variables and parameters to the object. This analysis is always for an isolated rotor, regardless of rc_type. If the rotorcraft type is coaxial, then modifications are made, as necessary.
        """

        if mode == "real":
            dtype = float
        elif mode == "complex":
            dtype = complex
        else:
            raise ValueError("Analysis mode must be either 'real' or 'complex'.")

        # Check to make sure that V_inf is greater than zero if flight mode is climb, descent, or forward
        V_inf = self.variables["V_inf"].value

        if (
            self.parameters["mode"] in ["climb", "descent", "forward"]
            and np.real(V_inf) > 0.0
        ):
            self.V_inf = V_inf
        elif self.parameters["mode"] == "hover" and V_inf == 0.0:
            self.V_inf = 0.0
        elif self.parameters["mode"] in ["axial_climb", "axial_descent"]:
            self.V_inf = V_inf  # For these modes, V_c = V_inf
        else:
            raise ValueError(
                f"Selected mode is '{self.parameters['mode']}', but the input values for this mode do not align with expected behavior. Verify that the selected mode matches the conditions specified for the set variables."
            )

        # Check to make sure that the rotor disk angle of attack is non-zero if fligth mode is climb, descent, or forward
        # alpha = self.variables["alpha"].value

        # if self.parameters["mode"] in ["climb", "descent", "forward"]:
        #     self.alpha = alpha
        # elif (
        #     self.parameters["mode"] in ["hover", "axial_climb", "axial_descent"]
        #     and alpha == 0.0
        # ):
        #     self.alpha = alpha
        # else:
        #     raise ValueError(
        #         f"Selected mode is {self.parameters['mode']}, but the input values for this mode do not align with expected behavior. Verify that the selected mode matches the conditions specified for the set variables."
        #     )

        # Adjust the value of alpha and compute other quantities if equilibrium conditions are specified
        self._equilibrium_check()

        # Compute the quantities of interest for hover
        self.C_T, self.lam_h, self.v_h = self._compute_hover_quantities()

        # Compute the advance ratio for the rotorcraft
        self.mu = (
            self.V_inf
            * np.cos(self.variables["alpha"].value)
            / (self.parameters["Omega"] * self.parameters["R"])
        )

        # Compute the vertical advance ratio for the rotorcraft
        self.mu_y = (
            self.V_inf
            * np.sin(self.variables["alpha"].value)
            / (self.parameters["Omega"] * self.parameters["R"])
        )

        if self.parameters["mode"] == "descent":
            self.mu_y *= -1.0

        # Compute the lambda (inflow)/power ratio using the hover quantities, advance ratio, and initial guess with the Newton-Raphson method for hover, climb, descent, or forward flight
        if self.parameters["mode"] in ["hover", "climb", "descent", "forward"]:
            # Set the initial guess for lam0 based on the mode
            if self.parameters["mode"] == "descent":
                lam0 = (
                    -self.V_inf
                    / self.parameters["Omega"]
                    / self.parameters["R"]
                    * np.tan(self.variables["alpha"].value)
                )
            else:
                lam0 = self.lam_h

            self.lam = self._compute_lam_newton(
                lam0=lam0,
                mu_x=self.mu,
                mu_y=self.mu_y,
                alpha=self.variables["alpha"].value,
                C_T=self.C_T,
            )

            # Compute the lambda/power ratio using the determined value for lambda
            self.LPR = self.lam / self.lam_h

            # Compute the power required to hover
            self.P_h = self.thrust * self.v_h

            # Compute the power for the system for induced and climb/propulsion
            P_tot = self.LPR * self.P_h

            # Compute the induced power for the rotor (i.e. isolate the induced power from the climb/propulsion power)
            self.P_i = self.lam_h * self.P_h / np.sqrt(self.mu**2 + self.lam**2)

            # # Compute the climb/propulsion power for the rotor
            self.P_cp = P_tot - self.P_i

        # Perform the alternative procedure for axial climb/descent
        else:
            # Set the climb velocity
            self.v_c = self.V_inf

            # Compute the ratio of climb velocity to induced hover velocity
            self.vc_vh = self.v_c / self.v_h

            # Compute the ratio of the induced velocity to the hover induced velocity
            self.vi_vh = self._eval_axial_climb_descent(v_c=self.v_c, v_h=self.v_h)

            # Compute the lambda/power ratio using the computed ratio for vi_vh and the climb velocity
            self.LPR = self.vc_vh + self.vi_vh

            # Compute the power required to hover
            self.P_h = self.thrust * self.v_h

            # Compute the induced power for the rotor (i.e. isolate the induced power from the climb/propulsion power)
            self.P_i = self.P_h * self.vi_vh

            # Compute the climb power for the rotor
            self.P_cp = self.vc_vh * self.P_h

        # Modify the induced and propulsion power by the appropriate factors if the rotorcraft type is coaxial
        if self.parameters["rc_type"] == "coaxial":
            # Compute the adjusted induced power
            P_i_coax = 2 * self.P_i * self.parameters["kappa_int"]
            self.P_i = P_i_coax

            # Compute the adjusted propulsion power
            P_cp_coax = 2 * self.P_cp * self.parameters["kappa_int"]
            self.P_cp = P_cp_coax

        # Compute the parasitic and profile power values with the known advance ratio
        self.P_0, self.P_p = self._compute_parasitic_profile_power()

        # Update the analyzed attribute
        self.analyzed = True

        # Assign the outputs to the outputs dictionary
        self.outputs = {}

        self.outputs["C_T"] = State(
            value=self.C_T, desc="Thrust coefficient for the rotor disk", source=self
        )

        self.outputs["lam_h"] = State(
            value=self.lam_h, desc="Inflow ratio in hover", source=self
        )

        self.outputs["v_h"] = State(
            value=self.v_h,
            desc="Induced velocity at the rotor disk in hover",
            source=self,
        )

        self.outputs["mu"] = State(
            value=self.mu,
            desc="Advance ratio of the rotorcraft for forward flight",
            source=self,
        )

        self.outputs["mu_y"] = State(
            value=self.mu_y,
            desc="Vertical advance ratio of the rotorcraft",
            source=self,
        )

        self.outputs["LPR"] = State(
            value=self.LPR,
            desc="Lambda (inflow)/power ratio for a single rotor in the system (i.e. lam/lam_h or P/P_h; this includes the induced power contributions and any power associated with climb/descent/forward flight). Note that if the system is coaxial, this output does not represent the ratio for the whole coaxial system.",
            source=self,
        )

        self.outputs["P_i"] = State(
            value=self.P_i,
            desc="Induced power for the rotor system (if rotorcraft type is coaxial, this induced power has already been modified to reflect the total induced power)",
            source=self,
        )

        self.outputs["P_cp"] = State(
            value=self.P_cp,
            desc="Power required for climb/propulsion of the rotorcraft.",
            source=self,
        )

        self.outputs["P_0"] = State(
            value=self.P_0,
            desc="Profile power coefficient for the rotor system (if rotorcraft type is coaxial, this has already been doubled to reflect the power for both rotor disks).",
            source=self,
        )

        self.outputs["P_p"] = State(
            value=self.P_p, desc="Parasitic power for the rotorcraft.", source=self
        )

        # Add the climb velocity ratio as an additional output if axial climb/descent is selected
        if self.parameters["mode"] in ["axial_climb", "axial_descent"]:
            self.outputs["vc_vh"] = State(
                value=self.vc_vh,
                desc="Climb velocity ratio (i.e. v_c/v_h) for the rotorcraft. This is only a provided output in the axial climb/descent operational modes.",
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
        self.C_Tb = self.outputs["C_T"].deriv
        self.lam_hb = self.outputs["lam_h"].deriv
        self.v_hb = self.outputs["v_h"].deriv
        self.mub = self.outputs["mu"].deriv
        self.mu_yb = self.outputs["mu_y"].deriv
        self.LPRb = self.outputs["LPR"].deriv
        self.P_ib = self.outputs["P_i"].deriv
        self.P_cpb = self.outputs["P_cp"].deriv
        self.P_0b = self.outputs["P_0"].deriv
        self.P_pb = self.outputs["P_p"].deriv

        if self.parameters["mode"] in ["axial_climb", "axial_descent"]:
            self.vc_vhb = self.outputs["vc_vh"].deriv

        # Compute the rotor tip velocity
        self.V_tip = self.parameters["Omega"] * self.parameters["R"]

        # PROFILE POWER
        # Multiply the bar value by 2 if the type is coaxial
        if self.parameters["rc_type"] == "coaxial":
            self.P_0b *= 2.0

        # Compute the derivative contribution to C_P0 from P_0
        self.C_P0b = (
            self.P_0b
            * self.parameters["rho"]
            * np.pi
            * self.parameters["R"] ** 2
            * self.V_tip**3
        )

        # Compute the derivative contribution to m_P0 from C_P0
        self.m_P0b = self.C_P0b * (1 + 4.65 * self.outputs["mu"].value ** 2)

        # Compute the derivative contribution to mu from C_P0
        self.mub += (
            self.C_P0b
            * (2 * 4.65 * self.outputs["mu"].value)
            * self.variables["m_P0"].value
        )

        # PARASITIC POWER
        # Compute the derivative contribution to C_Pp from P_p
        self.C_Ppb = (
            self.P_pb
            * self.parameters["rho"]
            * np.pi
            * self.parameters["R"] ** 2
            * self.V_tip**3
        )

        # Compute the derivative contribution to n_Pp from C_Pp
        self.n_Ppb = self.C_Ppb * self.outputs["mu"].value ** 3

        # Compute the derivative contribution to mu from C_Pp
        self.mub += (
            self.C_Ppb
            * self.variables["n_Pp"].value
            * 3
            * self.outputs["mu"].value ** 2
        )

        # INDUCED AND PROPULSION POWER
        self.V_infb = self.variables["V_inf"].deriv
        self.alphab = self.variables["alpha"].deriv

        self._induced_propulsion_power_adj()

        # ADVANCE RATIOS
        if self.parameters["mode"] in ["descent"]:
            self.mu_yb *= -1.0

        # Compute the contribution to V_inf from mu_y
        self.V_infb += self.mu_yb * np.sin(self.variables["alpha"].value) / self.V_tip

        # Compute the contribution to alpha from mu_y
        self.alphab += (
            self.mu_yb
            * self.variables["V_inf"].value
            * np.cos(self.variables["alpha"].value)
            / self.V_tip
        )

        # Compute the contribution to V_inf from mu
        self.V_infb += self.mub * np.cos(self.variables["alpha"].value) / self.V_tip

        self.alphab += (
            self.mub
            * -self.variables["V_inf"].value
            * np.sin(self.variables["alpha"].value)
            / self.V_tip
        )

        # HOVER QUANTITIES
        self._hover_quantities_adj()

        # EQUILIBRIUM CHECK CALCS
        if self.parameters["equilibrium"] and self.parameters["mode"] != "hover":

            # Compute the contribution to alpha from thrust if equilibrium conditions are specified
            self.alphab += (
                self.thrustb
                * self.weight
                * np.sin(self.variables["alpha"].value)
                / (np.cos(self.variables["alpha"].value) ** 2)
            )

            # Compute the contribution to D_f from alpha
            self.D_fb = self.alphab / (1 + (self.D_f / self.weight) ** 2) / self.weight

            # Compute the contribution to V_inf from D_f
            self.V_infb += (
                self.D_fb
                * self.parameters["rho"]
                * self.variables["V_inf"].value
                * self.f
            )

            # Compute the contribution to f from D_f
            self.fb = (
                self.D_fb
                * 0.5
                * self.parameters["rho"]
                * self.variables["V_inf"].value ** 2
            )

            # Compute the contribution to n_Pp from f
            self.n_Ppb += self.fb * 2.0 * np.pi * self.parameters["R"] ** 2

        # Set the derivative values for the variables
        self.variables["V_inf"].set_deriv_value(self.V_infb)
        self.variables["alpha"].set_deriv_value(self.alphab)
        self.variables["n_Pp"].set_deriv_value(self.n_Ppb)
        self.variables["m_P0"].set_deriv_value(self.m_P0b)

        # Update the attribute for the adjoint analysis
        self.adjoint_analyzed = True

        return

    def _compute_hover_quantities(self):
        """
        Comptues the hover quantities for the rotorcraft.

        Returns
        -------
        C_T : float
            Coefficient of thrust for the rotorcraft
        lam_h : float
            Inflow ratio in hover
        v_h : float
            Induced velocity in hover (m/s)
        """

        # Compute the coefficient of thrust, assuming lift = weight
        A = np.pi * self.parameters["R"] ** 2
        C_T = self.thrust / (
            self.parameters["rho"]
            * A
            * self.parameters["Omega"] ** 2
            * self.parameters["R"] ** 2
        )

        # Compute the inflow ratio in hover
        lam_h = np.sqrt(C_T / 2)

        # Compute the induced velocity in hover
        v_h = lam_h * self.parameters["Omega"] * self.parameters["R"]

        return C_T, lam_h, v_h

    def _compute_lam_newton(
        self, lam0, mu_x, mu_y, alpha, C_T, epsilon=1e-12, max_iter=100
    ):
        """
        Computes the inflow ratio for a rotorcraft using a Newton-Raphson procedure.

        Parameters
        ----------
        lam0 : float
            Initial guess for the inflow ratio
        mu_x : float
            Advance ratio in the horizontal direction (forward flight)
        mu_y : float
            Advance ratio in the vertical direction (climbing/descending flight)
        alpha : float
            Rotor disk angle of attack (rad)
        epsilon : float
            Tolerance for the Newton-Raphson solver, defaults to 1e-12
        max_iter : int
            Max number of iterations for the Newton-Raphson solver, defaults to 100 iterations

        Returns
        -------
        lam : float
            Converged value of the inflow ratio from solving the generalized inflow equation numerically
        """

        # Set the initial value for inflow to lam0
        lam = lam0

        # Loop for the Newton-Raphson procedure to solve for lambda
        for i in range(max_iter):

            # Compute the function/derivative values with the current data
            f_lam = self._eval_lam_func(lam, mu_x, mu_y, alpha, C_T)

            deriv_lam = self._eval_lam_deriv(lam, mu_x, C_T)

            # Compute the new value for lambda with the update
            lam_new = lam - f_lam / deriv_lam

            # Compute the error
            error = np.abs((lam_new - lam) / lam_new)

            # If the error is less than the tolerance or maximum number of iterations is reached, break; otherwise, continue with the procedure
            if error < epsilon:
                lam = lam_new
                return lam
            elif i == max_iter - 1:
                raise RuntimeError(
                    f"Maximum number of iterations reached, current error = {error:.8f}"
                )
            else:
                lam = lam_new
                continue

    def _eval_lam_func(self, lam, mu_x, mu_y, alpha, C_T):
        """
        Evaluates the residual equation for the generalized inflow equation.

        Parameters
        ----------
        lam : float
            Current value for the inflow ratio
        mu_x : float
            Advance ratio in the horizontal direction (forward flight)
        mu_y : float
            Advance ratio in the vertical direction (climbing/descending flight)
        alpha : float
            Rotor disk angle of attack (rad)
        C_T : float
            Coefficient of thrust for the rotor disk

        Returns
        -------
        f_lam : float
            Value of the residual equation using the provided inputs
        """

        # Evaluate the inflow ratio function (forward flight without climb)
        f_lam = lam - mu_x * np.tan(alpha) - C_T / (2 * np.sqrt(mu_x**2 + lam**2))

        # Use logic to determine whether mu_y contribution should be added (i.e. is the mode climb/descent)
        if self.parameters["mode"] in ["climb", "descent"]:
            f_lam -= mu_y

        return f_lam

    def _eval_lam_deriv(self, lam, mu_x, C_T):
        """
        Evaluates the derivative of the residual equation for the generalized inflow equation.

        Parameters
        ----------
        lam : float
            Current value for the inflow ratio
        mu_x : float
            Advance ratio in the horizontal direction (forward flight)
        C_T : float
            Coefficient of thrust for the rotor disk

        Returns
        -------
        deriv_lam : float
            Value of the derivative of the residual equation using the provided inputs
        """

        # Evaluate the derivative of the inflow ratio function
        deriv_lam = 1 + C_T / 2 * (mu_x**2 + lam**2) ** (-3 / 2) * lam

        return deriv_lam

    def _eval_axial_climb_descent(self, v_c, v_h):
        """
        Computes the induced velocity normalized by the induced velocity in hover for axial climb and descent. Value of v_c dictates which procedure is followed.

        Parameters
        ----------
        v_c : float
            Climb velocity for the rotorcraft (m/s)
        v_h : float
            Induced velocity in hover for the rotorcraft (m/s)

        Returns
        -------
        vi_vh : float
            Ratio of the induced velocity to the induced velocity in hover, where v_i is computed based on the value for v_c
        """

        # Check to make sure that the climb velocity is negative for axial descent
        if self.parameters["mode"] == "axial_descent" and v_c >= 0.0:
            raise ValueError(
                f"The selected mode is '{self.parameters['mode']}', and the value provided for the climb velocity is v_c = {v_c}. The climb velocity must be negative for axial descent."
            )

        # Evaluate the ratio for climb velocity to hover induced velocity
        C = v_c / (v_h)

        # Evaluate the rotor using normal equation when in normal working state for axial climb
        if self.parameters["mode"] == "axial_climb" and np.real(self.vc_vh) > 0.0:

            vi_vh = (-C / 2) + np.sqrt((C / 2) ** 2 + 1)
        # Evaluate the rotor using the modified equation when in windmill brake state
        elif self.parameters["mode"] == "axial_descent" and np.real(self.vc_vh) <= -2.0:

            vi_vh = (-C / 2) - np.sqrt((C / 2) ** 2 - 1)
        # If in the turbulent wake state or vortex ring state, use the continuous approximation, which is based on experimental data
        elif (
            self.parameters["mode"] == "axial_descent"
            and -2.0 < np.real(self.vc_vh) < 0.0
        ):
            kappa = self.parameters["kappa"]
            k1 = -1.125
            k2 = -1.372
            k3 = -1.718
            k4 = -0.655

            vi_vh = kappa + k1 * C + k2 * C**2 + k3 * C**3 + k4 * C**4
        else:
            raise ValueError(
                f"The selected mode '{self.parameters['mode']}' and the value for V_inf = V_c = {v_c} do not align. Verify that the selected mode and input conditions make physical sense."
            )

        return vi_vh

    def _equilibrium_check(self):
        """
        Sets the value for the thrust attribute depending on the Boolean value for the equilibrium parameter.

        Returns
        -------
        None, but internally sets the value for the thrust attribute
        """

        # Check to see if equilibrium conditions are specified
        if self.parameters["equilibrium"]:
            # Compute the equivalent flat plate area
            self.f = (
                self.variables["n_Pp"].value * 2.0 * np.pi * self.parameters["R"] ** 2
            )

            # Compute the parasitic drag
            self.D_f = (
                0.5
                * self.parameters["rho"]
                * self.variables["V_inf"].value ** 2
                * self.f
            )

            # Compute the rotor disk angle of attack using a force balance and known information
            alpha_eq = np.arctan2(np.real(self.D_f), self.weight)

            if isinstance(alpha_eq, np.ndarray):
                alpha_eq = alpha_eq[0]

            # Update the variable value for the rotor disk angle of attack
            self.set_var_values(variables={"alpha": alpha_eq})

            # Solve for the required thrust via the force equilibrium
            self.thrust = self.weight / np.cos(self.variables["alpha"].value)

        # If equilibrium conditions are not specified, just set thrust equal to weight
        else:
            self.thrust = self.weight

        return

    def _compute_parasitic_profile_power(self):
        """
        Computes the dimensional values of profile and parasitic power.

        Returns
        -------
        P_0 : float
            Dimensional value of the profile power for the rotorcraft (W)
        P_p : float
            Dimensional value of the parasitic power for the rotorcraft (W)
        """

        # Compute the profile power coefficient
        C_P0 = self.variables["m_P0"].value * (1 + 4.65 * self.mu**2)

        # Compute the dimensional profile power
        A = np.pi * self.parameters["R"] ** 2
        P_0 = (
            C_P0
            * self.parameters["rho"]
            * A
            * (self.parameters["Omega"] * self.parameters["R"]) ** 3
        )

        # If the rotorcraft type is coaxial, modify the profile power by multiplying by two
        if self.parameters["rc_type"] == "coaxial":
            P_0 *= 2.0

        # Compute the parasitic power coefficient
        C_Pp = self.variables["n_Pp"].value * self.mu**3

        # Compute the dimensional parasitic power
        P_p = (
            C_Pp
            * self.parameters["rho"]
            * A
            * (self.parameters["Omega"] * self.parameters["R"]) ** 3
        )

        return P_0, P_p

    def _induced_propulsion_power_adj(self):
        """
        Performs the adjoint implementation for the induced and propulsion power computations.
        """

        # If the rc_type is coaxial, multiply the bar values for P_cp and P_i by 2.0 and kappa_int
        if self.parameters["rc_type"] == "coaxial":
            self.P_cpb *= 2.0 * self.parameters["kappa_int"]
            self.P_ib *= 2.0 * self.parameters["kappa_int"]

        # Perform the adjoint procedure, where the operations are dependent on the operational mode chosen
        if self.parameters["mode"] in ["hover", "climb", "descent", "forward"]:
            # Compute the derivative contribution to P_tot from P_cp
            self.P_totb = self.P_cpb

            # Compute the derivative contribution to P_i from P_cp
            self.P_ib += -self.P_cpb

            # Compute the partial derivatives of the residual (implicit lambda) expression
            dRdCT = -1 / (2 * np.sqrt(self.mu**2 + self.lam**2))

            dRdmu = (
                -np.tan(self.variables["alpha"].value)
                + self.C_T / 2 * (self.mu**2 + self.lam**2) ** (-3 / 2) * self.mu
            )

            if self.parameters["mode"] in ["climb", "descent"]:
                dRdmuy = -1.0
            else:
                dRdmuy = 0.0

            dRdalpha = -self.mu / ((np.cos(self.variables["alpha"].value)) ** 2)

            # Compute the partial derivative of the residual expression wrt lambda
            dRdlam = (
                1 + self.C_T / 2 * (self.mu**2 + self.lam**2) ** (-3 / 2) * self.lam
            )

            # Compute the derivative contributions to the variables through the expression for P_i (note that there is explicit and implicit dependence here)
            dPidlam = (
                -self.lam_h
                * self.P_h
                * (self.mu**2 + self.lam**2) ** (-3 / 2)
                * self.lam
            )

            psi_Pi = (
                dRdlam ** (-1) * -dPidlam * self.P_ib
            )  # adjoint variables for P_i expression

            self.C_Tb += psi_Pi * dRdCT

            self.mub += (
                self.P_ib
                * -self.lam_h
                * self.P_h
                * (self.mu**2 + self.lam**2) ** (-3 / 2)
                * self.mu
                + psi_Pi * dRdmu
            )

            self.mu_yb += psi_Pi * dRdmuy

            self.alphab += psi_Pi * dRdalpha

            self.lam_hb += self.P_ib * self.P_h / np.sqrt(self.mu**2 + self.lam**2)

            self.P_hb = self.P_ib * self.lam_h / np.sqrt(self.mu**2 + self.lam**2)

            # Compute the derivative contribution to LPR from P_tot
            self.LPRb += self.P_totb * self.P_h

            # Compute the derivative contribution to P_h from P_tot
            self.P_hb += self.P_totb * self.LPR

            # Compute the derivative contribution to thrust from P_h
            self.thrustb = self.P_hb * self.v_h

            # Compute the derivative contribution to v_h from P_h
            self.v_hb += self.P_hb * self.thrust

            # Compute the derivative contributions to the variables dependent on lambda through the expression for LPR and the implicit function for lambda
            dLPRdlam = 1 / self.lam_h

            psi_LPR = (
                dRdlam ** (-1) * -dLPRdlam * self.LPRb
            )  # adjoint variables for LPR expression

            self.C_Tb += psi_LPR * dRdCT

            self.mub += psi_LPR * dRdmu

            self.mu_yb += psi_LPR * dRdmuy

            self.alphab += psi_LPR * dRdalpha

            self.lam_hb += self.LPRb * -self.lam / self.lam_h**2

        else:
            # Compute the derivative contribution to vc_vh from P_cp
            self.vc_vhb += self.P_cpb * self.P_h

            # Compute the derivative contribution to P_h from P_cp
            self.P_hb = self.P_cpb * self.vc_vh

            # Compute the derivative contribution to vi_vh from P_i
            self.vi_vhb = self.P_ib * self.P_h

            # Compute the derivative contribution to P_h from P_i
            self.P_hb += self.P_ib * self.vi_vh

            # Compute the derivative contribution to thrust from P_h
            self.thrustb = self.P_hb * self.v_h

            # Compute the derivative contribution to v_h from P_h
            self.v_hb += self.P_hb * self.thrust

            # Compute the derivative contribution to vc_vh from LPR
            self.vc_vhb += self.LPRb

            # Compute the derivative contribution to vi_vh from LPR
            self.vi_vhb += self.LPRb

            # Compute the derivative contribution to vc_vh from vi_vh
            if self.parameters["mode"] == "axial_climb" and self.vc_vh > 0.0:
                self.vc_vhb += self.vi_vhb * (
                    -0.5
                    + 0.5 * ((self.vc_vh / 2) ** 2 + 1) ** (-1 / 2) * self.vc_vh / 2
                )
            elif self.parameters["mode"] == "axial_descent" and self.vc_vh <= -2.0:
                self.vc_vhb += self.vi_vhb * (
                    -0.5
                    - 0.5 * ((self.vc_vh / 2) ** 2 - 1) ** (-1 / 2) * self.vc_vh / 2
                )
            elif self.parameters["mode"] == "axial_descent" and -2.0 < self.vc_vh < 0.0:
                k1 = -1.125
                k2 = -1.372
                k3 = -1.718
                k4 = -0.655

                self.vc_vhb += self.vi_vhb * (
                    k1
                    + k2 * 2 * self.vc_vh
                    + k3 * 3 * self.vc_vh**2
                    + k4 * 4 * self.vc_vh**3
                )
            else:
                raise ValueError(
                    f"The adjoint for these conditions is not defined. This should not occur, verify that inputs are correct for the operational mode and flight speed."
                )

            # Compute the derivative contribution to v_c from vc_vh
            self.v_cb = self.vc_vhb * 1 / self.v_h

            # Compute the derivative contribution to v_h from vc_vh
            self.v_hb += self.vc_vhb * -self.v_c / self.v_h**2

            # Compute the derivative contribution to V_inf from v_c
            self.V_infb += self.v_cb

        return

    def _hover_quantities_adj(self):
        """
        Performs the adjoint implementation for the hover quantity computations.
        """

        # Compute the contribution to lam_h from v_h
        self.lam_hb += self.v_hb * self.V_tip

        # Compute the contribution to C_T from lam_h
        self.C_Tb += self.lam_hb * 0.25 * (self.C_T / 2) ** (-1 / 2)

        # Compute the contribution to thrust from C_T
        self.thrustb += self.C_Tb / (
            self.parameters["rho"] * np.pi * self.parameters["R"] ** 2 * self.V_tip**2
        )

        return


if __name__ == "__main__":

    # Create the segment object
    seg = MomentumTheorySegment(seg_name="test", mode="descent", weight=20.0 * 9.81)

    # Set the variables
    seg.set_var_values(variables={"alpha": 0.0, "V_inf": 0.0})

    seg.analyze()

    outs = seg.get_output_values()
    for key in outs:
        print(f"{key} = {outs[key]}")
