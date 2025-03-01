"""
Coordinate systems
==================

In this module the base class for coordinate systems is defined along with several
most commonly used systems: Cartesian, cylindrical, and spherical systems.

If you want to add a custom coordinate system, subclass ``CurvilinearCoordinateSystem``
and define the scalar names and the Lamé coefficients related to it.

**Limitations:**

1. The systems share the same point of origin and are defined in relation to the same
   axes. The rotation and translation of coordinate systems are not supported.

2. The systems are orthogonal and right-handed.

3. The space is 3-dimensional.
"""

from __future__ import annotations

from abc import ABCMeta, abstractmethod

from sympy import Expr, Symbol, symbols, S, sin


class BaseCoordinateSystem(metaclass=ABCMeta):
    """
    Abstract base class for all coordinate systems.
    """

    _display_name: str
    _base_scalars: tuple[Symbol, Symbol, Symbol]
    _lame_coefficients: tuple[Expr, Expr, Expr]

    @classmethod
    @abstractmethod
    def scalar_names(cls) -> tuple[str, str, str]:
        """
        Names of the base scalars.
        """

    @abstractmethod
    def __init__(self, display_name: str) -> None:
        pass

    @property
    def display_name(self) -> str:
        """
        Display name of the coordinate system.
        """

        return self._display_name

    @property
    def base_scalars(self) -> tuple[Expr, Expr, Expr]:
        """
        A 3-tuple of symbols for base scalars.
        """

        return self._base_scalars

    @property
    def lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        """
        A 3-tuple of Lamé coefficients of the coordinate system.
        """

        return self._lame_coefficients

    @property
    def jacobian_determinant(self) -> Expr:
        """
        The determinant of the Jacobian matrix of transformation. Equivalently,
        the product of all the Lamé coefficients.
        """

        lame_ = self.lame_coefficients
        return lame_[0] * lame_[1] * lame_[2]

    def __str__(self) -> str:
        return self.display_name

    def __eq__(self, other: object) -> bool:
        return (
            type(other) is type(self)
            and self.lame_coefficients == other.lame_coefficients  # type: ignore[attr-defined]
        )

    def __getitem__(self, name: str) -> Expr:
        names = self.scalar_names()
        try:
            index = names.index(name)
        except ValueError as e:
            raise IndexError(f"Expect a scalar name from {names}, got {name!r}.") from e

        return self.base_scalars[index]


def generate_scalar_names(system: BaseCoordinateSystem) -> tuple[str, str, str]:
    """
    Subscript the scalar names with the given suffix.
    """

    suffix = system.display_name
    names = [f"{name}_{suffix}" for name in system.scalar_names()]
    return names[0], names[1], names[2]


class RectilinearCoordinateSystem(
    BaseCoordinateSystem
):  # pylint: disable=abstract-method
    """
    Abstract subclass for rectilinear coordinate systems.
    """


class CartesianCoordinateSystem(RectilinearCoordinateSystem):
    """
    Cartesian coordinate system.
    """

    _lame_coefficients = S.One, S.One, S.One

    @classmethod
    def scalar_names(cls) -> tuple[str, str, str]:
        return "x", "y", "z"

    def __init__(self, display_name: str) -> None:
        self._display_name = display_name

        x, y, z = symbols(
            generate_scalar_names(self),
            real=True,
        )

        self._base_scalars = x, y, z


class CurvilinearCoordinateSystem(
    BaseCoordinateSystem
):  # pylint: disable=abstract-method
    """
    Abstract subclass for curvilinear coordinate systems.
    """


class CylindricalCoordinateSystem(CurvilinearCoordinateSystem):
    """
    Cylindrical coordinate system.
    """

    @classmethod
    def scalar_names(cls) -> tuple[str, str, str]:
        return "rho", "phi", "z"

    def __init__(self, display_name: str) -> None:
        self._display_name = display_name

        rho_str, phi_str, z_str = generate_scalar_names(self)

        rho = symbols(rho_str, nonnegative=True)
        phi, z = symbols((phi_str, z_str), real=True)

        self._base_scalars = rho, phi, z
        self._lame_coefficients = S.One, rho, S.One


class SphericalCoordinateSystem(CurvilinearCoordinateSystem):
    """
    Spherical coordinate system.
    """

    @classmethod
    def scalar_names(cls) -> tuple[str, str, str]:
        return "r", "theta", "phi"

    def __init__(self, display_name: str) -> None:
        self._display_name = display_name

        r_str, theta_str, phi_str = generate_scalar_names(self)

        r = symbols(r_str, nonnegative=True)
        theta, phi = symbols((theta_str, phi_str), real=True)

        self._base_scalars = r, theta, phi
        self._lame_coefficients = S.One, r, r * sin(theta)


__all__ = [
    "BaseCoordinateSystem",
    "RectilinearCoordinateSystem",
    "CartesianCoordinateSystem",
    "CurvilinearCoordinateSystem",
    "CylindricalCoordinateSystem",
    "SphericalCoordinateSystem",
    "generate_scalar_names",
]
