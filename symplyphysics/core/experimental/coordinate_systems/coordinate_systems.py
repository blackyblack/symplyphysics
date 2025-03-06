from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import Any, Optional, Sequence, Self
from sympy import Basic, Expr, S, sin

from symplyphysics import Symbol, angle_type, units
from symplyphysics.core.dimensions import dimensionless, Dimension, assert_equivalent_dimension
from ..vectors import VectorExpr, VectorSymbol


def _check_base_scalars(
    base_scalars: Sequence[Symbol],
    expected_dimensions: Sequence[Dimension],
) -> None:
    for base_scalar, expected_dimension in zip(base_scalars, expected_dimensions):
        assert_equivalent_dimension(
            base_scalar.dimension,
            "base_scalar",
            "_check_base_scalars",
            expected_dimension,
        )


def _check_base_vectors(base_vectors: Sequence[VectorSymbol]) -> None:
    for base_vector in base_vectors:
        assert_equivalent_dimension(
            base_vector.dimension,
            "base_vector",
            "_check_base_vectors",
            dimensionless,
        )

        norm = base_vector.norm
        if norm != 1:
            raise ValueError(f"Expected {base_vector} to be unit, but its norm is {norm}.")

    # Since all base_vectors are symbols and not expressions, we can check their linear
    # independence by simple equality
    base_vector_set = set(base_vectors)
    if len(base_vectors) != len(base_vector_set):
        raise ValueError(f"Expected 3 base vectors, got {base_vector_set}.")


class BaseCoordinateSystem(Basic, metaclass=ABCMeta):  # type: ignore[misc]
    """
    Base class for coordinate systems.

    Note that we only deal with **orthogonal** coordinates, in which the coordinate surfaces meet
    at right angles at all points in space. Only in such coordinates the metric tensor is diagonal,
    swhich simplifies the calculations. The opposite is called **skew** coordinates.
    """

    def __new__(
        cls,
        *,
        base_scalars: Optional[Sequence[Symbol]] = None,  # pylint: disable=unused-argument
        base_vectors: Optional[Sequence[VectorSymbol]] = None  # pylint: disable=unused-argument
    ) -> Self:
        return Basic.__new__(cls)  # type: ignore[no-any-return]

    def __init__(
        self,
        *,
        base_scalars: Optional[Sequence[Symbol]] = None,
        base_vectors: Optional[Sequence[VectorSymbol]] = None,
    ) -> None:
        if base_scalars:
            s1, s2, s3 = base_scalars
            _check_base_scalars(base_scalars, self._base_scalar_dimensions())
        else:
            s1, s2, s3 = self._generate_base_scalars()

        base_scalars = s1, s2, s3

        if base_vectors:
            v1, v2, v3 = base_vectors
            _check_base_vectors(base_vectors)
        else:
            v1, v2, v3 = self._generate_base_vectors()

        base_vectors = v1, v2, v3

        self._args = base_scalars, base_vectors

    @staticmethod
    @abstractmethod
    def _base_scalar_dimensions() -> tuple[Dimension, Dimension, Dimension]:
        pass

    @staticmethod
    @abstractmethod
    def _generate_base_scalars() -> tuple[Symbol, Symbol, Symbol]:
        pass

    @staticmethod
    @abstractmethod
    def _generate_base_vectors() -> tuple[VectorSymbol, VectorSymbol, VectorSymbol]:
        pass

    @property
    def base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return self.args[0]  # type: ignore[no-any-return]

    @property
    def base_vectors(self) -> tuple[VectorSymbol, VectorSymbol, VectorSymbol]:
        return self.args[1]  # type: ignore[no-any-return]

    @property
    @abstractmethod
    def lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        pass

    @property
    def jacobian(self) -> Expr:
        h1, h2, h3 = self.lame_coefficients
        return h1 * h2 * h3

    def _eval_subs(self, *args: Any, **kwargs: Any) -> BaseCoordinateSystem:
        # NOTE that the coordinate system object itself can't be substituted.

        base_scalars = [scalar.subs(*args, **kwargs) for scalar in self.base_scalars]
        base_vectors = [vector.subs(*args, **kwargs) for vector in self.base_vectors]
        return type(self)(base_scalars=base_scalars, base_vectors=base_vectors)


class RectilinearCoordinateSystem(BaseCoordinateSystem):  # pylint: disable=abstract-method
    """
    A **rectilinear** coordinate system is based on the intersection of straight lines.
    """


class CartesianCoordinateSystem(RectilinearCoordinateSystem):

    @staticmethod
    def _base_scalar_dimensions() -> tuple[Dimension, Dimension, Dimension]:
        return units.length, units.length, units.length

    @staticmethod
    def _generate_base_scalars() -> tuple[Symbol, Symbol, Symbol]:
        x = Symbol("x", units.length, real=True)
        y = Symbol("y", units.length, real=True)
        z = Symbol("z", units.length, real=True)

        return x, y, z

    @staticmethod
    def _generate_base_vectors() -> tuple[VectorSymbol, VectorSymbol, VectorSymbol]:
        i = VectorSymbol("i", norm=1, display_latex="\\hat{\\mathbf{\\imath}}")
        j = VectorSymbol("j", norm=1, display_latex="\\hat{\\mathbf{\\jmath}}")
        k = VectorSymbol("k", norm=1, display_latex="\\hat{\\mathbf{k}}")

        return i, j, k

    @property
    def x(self) -> Symbol:
        return self.base_scalars[0]

    @property
    def y(self) -> Symbol:
        return self.base_scalars[1]

    @property
    def z(self) -> Symbol:
        return self.base_scalars[2]

    @property
    def i(self) -> VectorSymbol:
        return self.base_vectors[0]

    @property
    def j(self) -> VectorSymbol:
        return self.base_vectors[1]

    @property
    def k(self) -> VectorExpr:
        return self.base_vectors[2]

    @property
    def lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        return S.One, S.One, S.One


class CurvilinearCoordinateSystem(BaseCoordinateSystem):  # pylint: disable=abstract-method
    """
    A **curvilinear** coordinate system (in a narrow sense) is based on the intersection of curved
    lines, or a combination of curved and straight lines.
    """


class CylindricalCoordinateSystem(CurvilinearCoordinateSystem):

    @staticmethod
    def _base_scalar_dimensions() -> tuple[Dimension, Dimension, Dimension]:
        return units.length, angle_type, units.length

    @staticmethod
    def _generate_base_scalars() -> tuple[Symbol, Symbol, Symbol]:
        rho = Symbol("rho", units.length, display_latex="\\rho", nonnegative=True)
        phi = Symbol("phi", angle_type, display_latex="\\varphi", real=True)
        z = Symbol("z", units.length, real=True)

        return rho, phi, z

    @staticmethod
    def _generate_base_vectors() -> tuple[VectorSymbol, VectorSymbol, VectorSymbol]:
        e_rho = VectorSymbol("e_rho", norm=1, display_latex="\\hat{\\mathbf{e}}_{\\rho}")
        e_phi = VectorSymbol("e_phi", norm=1, display_latex="\\hat{\\mathbf{e}}_{\\varphi}")
        e_z = VectorSymbol("e_z", norm=1, display_latex="\\hat{\\mathbf{e}}_z")

        return e_rho, e_phi, e_z

    @property
    def rho(self) -> Symbol:
        return self.base_scalars[0]

    @property
    def phi(self) -> Symbol:
        return self.base_scalars[1]

    @property
    def z(self) -> Symbol:
        return self.base_scalars[2]

    @property
    def lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        return S.One, self.rho, S.One


class SphericalCoordinateSystem(CurvilinearCoordinateSystem):

    @staticmethod
    def _base_scalar_dimensions() -> tuple[Dimension, Dimension, Dimension]:
        return units.length, angle_type, angle_type

    @staticmethod
    def _generate_base_scalars() -> tuple[Symbol, Symbol, Symbol]:
        r = Symbol("r", units.length, nonnegative=True)
        theta = Symbol("theta", angle_type, display_latex="\\theta", nonnegative=True)
        phi = Symbol("phi", angle_type, display_latex="\\varphi", real=True)

        return r, theta, phi

    @staticmethod
    def _generate_base_vectors() -> tuple[VectorSymbol, VectorSymbol, VectorSymbol]:
        e_r = VectorSymbol("e_r", norm=1, display_latex="\\hat{\\mathbf{e}}_r")
        e_theta = VectorSymbol("e_theta", norm=1, display_latex="\\hat{\\mathbf{e}}_{\\theta}")
        e_phi = VectorSymbol("e_phi", norm=1, display_latex="\\hat{\\mathbf{e}}_{\\varphi}")

        return e_r, e_theta, e_phi

    @property
    def r(self) -> Symbol:
        return self.base_scalars[0]

    @property
    def theta(self) -> Symbol:
        return self.base_scalars[1]

    @property
    def phi(self) -> Symbol:
        return self.base_scalars[2]

    @property
    def lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        r = self.r
        theta = self.theta

        h_r = S.One
        h_theta = r
        h_phi = r * sin(theta)

        return h_r, h_theta, h_phi


__all__ = [
    "BaseCoordinateSystem",
    "RectilinearCoordinateSystem",
    "CartesianCoordinateSystem",
    "CurvilinearCoordinateSystem",
    "CylindricalCoordinateSystem",
    "SphericalCoordinateSystem",
]
