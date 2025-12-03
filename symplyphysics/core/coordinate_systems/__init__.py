from .base_system import BaseCoordinateSystem
from .cartesian_system import CartesianCoordinateSystem
from .cylindrical_system import CylindricalCoordinateSystem
from .spherical_system import SphericalCoordinateSystem

from .point import AppliedPoint
from .scalar import CoordinateScalar
from .vector import CoordinateVector, QuantityCoordinateVector, combine_coordinate_vectors

CARTESIAN = CartesianCoordinateSystem()
CYLINDRICAL = CylindricalCoordinateSystem()
SPHERICAL = SphericalCoordinateSystem()

__all__ = [
    # Most common coordinate system classes
    "BaseCoordinateSystem",
    "CartesianCoordinateSystem",
    "CylindricalCoordinateSystem",
    "SphericalCoordinateSystem",

    # Global constants
    "CARTESIAN",
    "CYLINDRICAL",
    "SPHERICAL",

    # Miscellaneous
    "AppliedPoint",
    "CoordinateScalar",
    "CoordinateVector",
    "QuantityCoordinateVector",
    "combine_coordinate_vectors",
]
