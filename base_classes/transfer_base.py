import numpy as np
import warnings


# TODO: need to clarify the name here
class TransferBase:
    """
    FIXME:
    """

    def __init__(self, obj_name: str, analysis_obj, **kwargs):

        # Store the name for the object
        self.obj_name = obj_name

        # Store the analysis object
        self.analysis_obj = analysis_obj

        # Set the parameters for the object using the kwargs or default parameters
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

        # Perform the check for the required methods for the object
        self.has_methods = self._check_methods()

        # Set the attribute that determines whether the adjoint has been initialized
        self.adjoint_initialized = False

        # Set the attribute that determines whether the data has been transferred from the analysis object
        self.transferred = False

    def _check_methods(self) -> bool:
        """
        FIXME:
        """

        # Check for the _transfer_analysis method
        if not (
            hasattr(self.__class__, "_transfer_analysis")
            and callable(getattr(self.__class__, "_transfer_analysis"))
        ):
            raise AttributeError(
                f"The class {self.__class__.__name__} does not have a method named '_transfer_analysis', which is a required method. \n"
            )

        has_methods = True

        return has_methods

    def _transfer_data(self):
        """
        FIXME:
        """

        # Check to see if the analysis object has been analyzed yet; if it has not, call the analyze method
        if not hasattr(self.analysis_obj, "stack"):
            self.analysis_obj.analyze()
        else:

            for obj in self.analysis_obj.stack:
                # Check to see if each object in the stack has been analyzed
                if not hasattr(obj, "analyzed") or not obj.analyzed:
                    # If one of the above conditions fails, perform the analysis procedure again for the whole stack (to make sure changes in variables propagate through)
                    self.analysis_obj.analyze()

        # if not hasattr(self.analysis_obj, "analyzed") or not self.analysis_obj.analyzed:
        #     self.analysis_obj.analyze()

        # Transfer the data from the analysis object's outputs to the current object
        self.obj_data = self.analysis_obj.outputs

        # Update the object that specifies whether or not data has been transferred
        self.transferred = True

        return

    def _parameter_initialization(self, kwargs):
        """
        FIXME:
        """

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

        return