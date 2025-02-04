import numpy as np
from flume.base_classes.analysis import Analysis
from examples.momentum_theory.momentum_model import MomentumTheorySegment
from examples.momentum_theory.profile_parasitic_power import ProfileParasiticPower

if __name__ == "__main__":
    # ---------- ISOLATED ANALYSIS MODEL ----------
    isolated = False
    if isolated:
        # Create the MomentumTheorySegment object by providing a name, operational mode, and specifying any sub-analyses
        hover = MomentumTheorySegment(
            seg_name="hover_model", mode="hover", sub_analysis=None
        )

        # The description of kwargs can be seen during inititialization (redefine above with desired input parameters)
        hover = MomentumTheorySegment(
            seg_name="hover_model",
            mode="hover",
            sub_analysis=None,
            Omega=30.0,
            R=2.0,
            weight=20.0,
            rc_type="coaxial",
        )

        # NOTE: comment out to avoid error
        # # If an incorrect kwarg is provided, error will be raised
        # hover = MomentumTheorySegment(
        #     seg_name="hover_model",
        #     mode="hover",
        #     sub_analysis=None,
        #     rpm=30.0,  # This is NOT a valid input parameter
        #     R=2.0,
        #     weight=20.0,
        #     rc_type="coaxial",
        # )

        # To display the variable information, call the metadata function (default returns info for all variables, but you can also specify a subset). Metadata is assigned to a dictionary, but can also be displayed by setting display_info to True
        # TODO: add function to display meta data dictionary
        var_meta = hover.get_var_metadata(display_info=True)

        # To extract and display variable/parameter values, use the get_var/parameter_values method
        param_values = hover.get_parameter_values(display_data=True)
        var_values = hover.get_var_values(display_data=True)

        # To change variable values form the default, use the set_var_values method; this will also raise an error if variable does not exist
        alpha = 0.0
        V_inf = 0.0
        m_P0 = 10.0

        # NOTE: comment out to avoid error
        # hover.set_var_values(variables={"AOA": alpha, "V_inf": V_inf})

        hover.set_var_values(variables={"alpha": alpha, "V_inf": V_inf, "m_P0": m_P0})

        # NOTE: comment out to avoid error
        # Output metadata/values can be extracted AFTER performing the analysis
        # hover.get_output_metadata()

        # To analyze the system, call the anlyze method
        hover.analyze()

        # Now, extract output metadata/values using similar methods as with parameter/variables
        out_meta = hover.get_output_metadata(display_info=True)
        out_vars = hover.get_output_values(display_data=True)

        # To declare design variables, which is used for housekeeping purposes and adjoint tests, use the declare_design_vars method
        hover.declare_design_vars(variables=["m_P0", "n_Pp"])

        # To test adjoint procedures for an analysis class, use the isolated adjoint procedure (to see debugging statements, use debug_print=True)
        hover.test_isolated_adjoint(debug_print=False)

    # ---------- COMBINED ANALYSIS MODEL ----------
    combined = True
    if combined:
        # Create the profile parasitic power model (use default variable values)
        ppp_model = ProfileParasiticPower(obj_name="ppp_model", R=2.0)

        # Create the axial climb model with the PPP model specified as a sub-analysis
        ax_climb = MomentumTheorySegment(
            seg_name="axclimb_model",
            mode="axial_climb",
            sub_analysis=ppp_model,
            Omega=30.0,
            R=2.0,
            weight=20.0,
            rc_type="coaxial",
        )

        ax_climb.set_var_values(variables={"V_inf": 2.0})

        # To analyze the entire system, simply call analyze on the top-level model. Internally, a stack is created and each analysis object is analyzed in sequence (executes from lowest-level sub-analysis up to the top)
        ax_climb.analyze()

        # Since ppp_model has outputs that are named the same as variables for hover model, they are connected during the analysis procedure. This is where a State's 'source' info is used
        connect_str = f"\n Connections map for object named '{ax_climb.obj_name}':"
        print(connect_str)
        print("-" * len(connect_str))
        for key in ax_climb.connects:
            print(
                f"  {key} : {ax_climb.connects[key]}"
            )  # for demo purposes only, this is used internally

        # Adjoint test can also be called for
        ppp_model.declare_design_vars(variables=["f", "C_d0", "sigma"])
        ax_climb.declare_design_vars(variables=["V_inf"])

        ax_climb.test_combined_adjoint(debug_print=False)
