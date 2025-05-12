from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import Any, Optional, Sequence, Self
from sympy import Basic, Expr, S, sin
from sympy.combinatorics.permutations import Permutation

from symplyphysics import Symbol, angle_type, units
from symplyphysics.core.dimensions import dimensionless, Dimension, assert_equivalent_dimension
from ..vectors import VectorSymbol, VectorFunction, AppliedVectorFunction, VectorNorm, VectorDot, VectorCross


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


def _check_base_vectors(base_vectors: Sequence[VectorSymbol | VectorFunction]) -> None:
    for base_vector in base_vectors:
        assert_equivalent_dimension(
            base_vector.dimension,
            "base_vector",
            "_check_base_vectors",
            dimensionless,
        )

        # TODO: add a check to see if `base_vector` has unit length

    # Since all base_vectors are symbols and not expressions, we can check their linear
    # independence by simple equality
    base_vector_set = set(base_vectors)
    if len(base_vectors) != len(base_vector_set):
        raise ValueError(f"Expected 3 base vectors, got {base_vector_set}.")


class BaseCoordinateSystem(Basic, metaclass=ABCMeta):
    """
    Base class for coordinate systems.

    Note that we only deal with **orthogonal** coordinates, in which the coordinate surfaces meet
    at right angles at all points in space. Only in such coordinates the metric tensor is diagonal,
    which simplifies the calculations. The opposite is called **skew** coordinates.
    """

    def __new__(cls,
        *,
        base_scalars: Optional[Sequence[Symbol]] = None,
        base_vectors: Optional[Sequence[VectorSymbol | VectorFunction]] = None) -> Self:
        return Basic.__new__(cls)

    def __init__(
        self,
        *,
        base_scalars: Optional[Sequence[Symbol]] = None,
        base_vectors: Optional[Sequence[VectorSymbol | VectorFunction]] = None,
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
    def _generate_base_vectors() -> tuple[
        VectorSymbol | VectorFunction,
        VectorSymbol | VectorFunction,
        VectorSymbol | VectorFunction,
    ]:
        pass

    @property
    def base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return self.args[0]

    def base_vectors(
        self,
        *args: Any,
    ) -> tuple[
            VectorSymbol | AppliedVectorFunction,
            VectorSymbol | AppliedVectorFunction,
            VectorSymbol | AppliedVectorFunction,
    ]:

        def prepare(
            base_vector: VectorSymbol | VectorFunction,) -> VectorSymbol | AppliedVectorFunction:
            if callable(base_vector):
                return base_vector(*args)

            return base_vector

        v1, v2, v3 = self.args[1]

        return prepare(v1), prepare(v2), prepare(v3)

    @property
    @abstractmethod
    def lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        pass

    @property
    def jacobian(self) -> Expr:
        h1, h2, h3 = self.lame_coefficients
        return h1 * h2 * h3

    def expr_simplify(self, expr: Expr, *args: Expr) -> Expr:
        """
        Simplifies the norm, dot and cross products of the base vectors of `self` in `expr`.

        `*args` is used to instantiate base vectors that are not coordinate-independent.
        """

        vs = self.base_vectors(*args)

        expr = expr.simplify()

        for i, v_i in enumerate(vs):
            expr = expr.subs(VectorNorm(v_i), 1)

            for j, v_j in enumerate(vs):
                if i == j:
                    dot = 1
                    cross = 0
                else:
                    dot = 0

                    ijk = {0, 1, 2}
                    ijk.remove(i)
                    ijk.remove(j)
                    (k,) = ijk

                    sign = Permutation((i, j, k)).signature()

                    cross = sign * vs[k]

                expr = expr.subs({
                    VectorDot(v_i, v_j, evaluate=False): dot,
                    VectorCross(v_i, v_j, evaluate=False): cross,
                })

        return expr


class CartesianCoordinateSystem(BaseCoordinateSystem):

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
        i = VectorSymbol("i", display_latex="\\hat{\\mathbf{\\imath}}")
        j = VectorSymbol("j", display_latex="\\hat{\\mathbf{\\jmath}}")
        k = VectorSymbol("k", display_latex="\\hat{\\mathbf{k}}")

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
        return self.base_vectors()[0]

    @property
    def j(self) -> VectorSymbol:
        return self.base_vectors()[1]

    @property
    def k(self) -> VectorSymbol:
        return self.base_vectors()[2]

    @property
    def lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        return S.One, S.One, S.One


class CylindricalCoordinateSystem(BaseCoordinateSystem):

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
    def _generate_base_vectors() -> tuple[VectorFunction, VectorFunction, VectorSymbol]:
        e_rho = VectorFunction("e_rho", display_latex="\\hat{\\mathbf{e}}_{\\rho}", nargs=1)
        e_phi = VectorFunction("e_phi", display_latex="\\hat{\\mathbf{e}}_{\\varphi}", nargs=1)
        e_z = VectorSymbol("e_z", display_latex="\\hat{\\mathbf{e}}_z")

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


class SphericalCoordinateSystem(BaseCoordinateSystem):

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
    def _generate_base_vectors() -> tuple[VectorFunction, VectorFunction, VectorFunction]:
        # TODO: replace with BaseVector later
        e_r = VectorFunction("e_r", display_latex="\\hat{\\mathbf{e}}_r", nargs=1)
        e_theta = VectorFunction("e_theta", display_latex="\\hat{\\mathbf{e}}_{\\theta}", nargs=1)
        e_phi = VectorFunction("e_phi", display_latex="\\hat{\\mathbf{e}}_{\\varphi}", nargs=1)

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
    "CartesianCoordinateSystem",
    "CylindricalCoordinateSystem",
    "SphericalCoordinateSystem",
]
