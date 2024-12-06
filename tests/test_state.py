import numpy as np
from base_classes.state import State
import unittest


class TestStateConstruction(unittest.TestCase):
    """
    Tests the construction of a State object from the framework.
    """

    def setUp(self):
        """
        Sets up the sample State object that is used for the subsequent tests.
        """

        self.scalar_state = State(value=2.0, desc="Scalar state")

        val = np.array([1.0, 2.0, 3.0])
        self.vec_state = State(value=val, desc="Vector state")

        return

    def test_scalar_value(self):
        """
        Checks the value of the scalar state
        """

        self.assertEqual(2.0, self.scalar_state.value, "Values do not match")

        return

    def test_scalar_desc(self):
        """
        Checks the description of the scalar state
        """

        test_str = "Scalar state"
        self.assertEqual(test_str, self.scalar_state.desc, "Descriptions do not match")

        return

    def test_scalar_type(self):
        """
        Checks the type of the scalar state value
        """

        self.assertEqual(self.scalar_state.data_type, float, "Types do not match")

        return

    def test_vec_value(self):
        """
        Checks the value of the vector state, which uses numpy
        """

        # self.assertEqual(np.array([1.0, 2.0, 3.0]), self.vec_state.value, "Values do not match")
        np.testing.assert_array_equal(
            np.array([1.0, 2.0, 3.0]), self.vec_state.value, "Values do not match"
        )

        return

    def test_vec_desc(self):
        """
        Checks the description of the vector state
        """

        test_str = "Vector state"
        self.assertEqual(test_str, self.vec_state.desc, "Descriptions do not match")

        return

    def test_vec_type(self):
        """
        Checks the type of the vector state
        """

        self.assertEqual(self.vec_state.data_type, np.ndarray, "Types do not match")

        return

    def test_vec_shape(self):
        """
        Checks the shape of the vector state
        """

        self.assertEqual(self.vec_state.shape, (3,), "Shapes do not match")


if __name__ == "__main__":
    # Run the unittests
    unittest.main()
