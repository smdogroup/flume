import numpy as np


class State:
    """
    This is a class to use for establishing a structure for variables and outputs of analysis classes.
    """

    def __init__(self, value, desc: str, deriv=None, source=None):
        """
        Base class that is used to wrap variables and outputs within the Flume framework.

        Parameters
        ----------
        value : float or np.ndarray
            Numeric value for the State
        desc : str
            String that describes the State
        deriv : float or np.ndarray
            Value for the State's derivative. This is None by default, and the user does not need to provide a numeric value when constructing the object
        source : class instance
            Object where the State data is sourced from. When creating states, this should be set to "self", otherwise the framework will raise an error

        Attributes
        ----------
        data_type : type
            The data type for the value of the State
        shape : tuple
            This is only created if the value for the State is a NumPy array, and it corresponds to the shape for the value
        """
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
        else:
            raise RuntimeError(
                "Argument 'source' was not set! Make sure that this is set to 'self' during default State construction within the __init__ method."
            )

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
