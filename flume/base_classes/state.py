import numpy as np


class State:
    """
    This is a class to use for establishing a structure for variables and outputs of analysis classes.

    Storing source information here will establish the flow of data between objects. Then, to access variable/output values one would do var[name].value; could do the same for source, deriv, etc.
    """

    def __init__(self, value, desc: str, deriv=None, source=None):
        # Set the variable value
        self.value = value

        # Set the variable data type
        self.data_type = type(self.value)

        # Set the variable shape if the value is a numpy array
        if isinstance(self.value, np.ndarray):
            self.shape = np.shape(self.value)

        # Set the description for the variable
        self.desc = desc

        # Set the derivative value
        # self.deriv = deriv
        if deriv is not None:
            self.deriv = deriv
        else:
            self.deriv = np.zeros_like(self.value)

        # Set the source for the variable
        if source is not None:
            self.source = source

    def set_deriv_value(self, deriv_val):
        """
        Helper method to set the derivative value for a State object.

        Parameters
        ----------
        deriv_val : float or numpy.ndarray
            The value fo the derivative for the State object.

        Returns
        -------
        None
        """

        self.deriv = deriv_val

        return
