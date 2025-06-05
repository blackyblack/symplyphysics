from __future__ import annotations

from typing import Optional, Any, Iterable, Self

from sympy import (ImmutableMatrix, Expr, sqrt, Symbol as SymSymbol, Basic, Derivative as
    SymDerivative, S)
from sympy.matrices.dense import DenseMatrix
from sympy.physics.units import Dimension
from sympy.physics.units.systems.si import dimsys_SI

from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.dimensions.collect_quantity import collect_quantity_factor_and_dimension
from symplyphysics.core.dimensions.miscellaneous import is_any_dimension, dimensionless
from symplyphysics.core.errors import UnitsError

from ..miscellaneous import sympify_expr
from ..vectors import VectorExpr, is_vector_expr, into_terms, split_factor
from .base_system import BaseCoordinateSystem
from .cartesian_system import CartesianCoordinateSystem
from .point import AppliedPoint, check_point_with_system, GLOBAL_POINT


class CoordinateVector(VectorExpr):
    """
    A `CoordinateVector` represents a vector in a certain coordinate _system_ at a particular
    _point_ in space and is defined as a linear combination of the system's base vectors. The
    coefficients of the linear combination are the *components* of the vector.

    Note that if a vector is defined in a Cartesian coordinate system, it will be equal to another
    component-based vector as long as their components match, regardless of their points of
    application. This is unlike vectors defined in non-Cartesian systems, which are not
    transferrable and therefore must have the same point of application in order to be compared.
    """

    _components: ImmutableMatrix
    _system: BaseCoordinateSystem
    _point: AppliedPoint | SymSymbol

    is_number = False
    is_atomic_vector = True

    @property
    def components(self) -> ImmutableMatrix:
        return self._components

    @property
    def system(self) -> BaseCoordinateSystem:
        return self._system

    @property
    def point(self) -> AppliedPoint | SymSymbol:
        return self._point

    @property
    def args(self) -> tuple[Basic, ...]:
        return self.components, self.system, self.point

    @property
    def free_symbols(self) -> set[Expr]:
        return self.components.free_symbols  # type: ignore[no-any-return]

    def _hashable_content(self) -> tuple[Basic, ...]:
        if isinstance(self.system, CartesianCoordinateSystem):
            # Cartesian-defined vectors are independent of their point of application.
            point = GLOBAL_POINT
        else:
            point = self.point

        return self.components, self.system, point

    def __new__(
        cls,
        components: Iterable[Any],
        system: BaseCoordinateSystem,
        point: Optional[AppliedPoint | SymSymbol] = None,
    ) -> Expr:
        if isinstance(components, DenseMatrix) and 1 not in components.shape:
            rows, cols = components.shape
            message = f"Expected a row or column vector, got a {rows} by {cols} matrix"
            raise ValueError(message)

        if not isinstance(components, DenseMatrix):
            components = ImmutableMatrix([sympify_expr(component) for component in components])
        else:
            components = ImmutableMatrix(components).reshape(3, 1)

        if all(component == 0 for component in components):
            return S.Zero

        obj = super().__new__(cls)  # pylint: disable=no-value-for-parameter

        point = check_point_with_system(system, point)

        obj._components = components
        obj._system = system
        obj._point = point

        return obj

    def __iter__(self) -> Iterable[Expr]:
        return iter(self.components)

    def __str__(self) -> str:
        name = type(self.system).__name__.removesuffix("CoordinateSystem")

        inner = ", ".join(str(component) for component in self.components)

        return f"{name}Vector({inner})"

    def _eval_vector_norm(self) -> Expr:
        return sqrt(sum(elem**2 for elem in self.components))

    @classmethod
    def _eval_vector_dot(cls, lhs: VectorExpr, rhs: VectorExpr) -> Optional[Expr]:
        if not (isinstance(lhs, CoordinateVector) and isinstance(rhs, CoordinateVector)):
            return None

        if type(lhs.system) is not type(rhs.system):
            return None

        return lhs.components.dot(rhs.components)

    @classmethod
    def _eval_vector_cross(cls, lhs: VectorExpr, rhs: VectorExpr) -> Optional[VectorExpr]:
        if not (isinstance(lhs, CoordinateVector) and isinstance(rhs, CoordinateVector)):
            return None

        if type(lhs.system) is not type(rhs.system) or lhs.point != rhs.point:
            return None

        components = lhs.components.cross(rhs.components)

        return CoordinateVector(components, lhs.system, lhs.point)

    def _eval_derivative(self, symbol: Expr) -> VectorExpr:
        # Needed for `sympy.solve` to work
        if is_vector_expr(symbol):
            return SymDerivative(self, symbol, evaluate=False)  # type: ignore[no-any-return]

        diff_components = ImmutableMatrix(
            [component.diff(symbol) for component in self.components],).T

        if isinstance(self.system, CartesianCoordinateSystem):
            diff_base_vectors = ImmutableMatrix.zeros(1, 3)
        else:
            applied_base_scalars = [f(symbol) for f in self.system.base_scalar_functions]
            diff_matrix = self.system.diff_base_vector_matrix(symbol, applied_base_scalars)
            diff_base_vectors = self.components.T * diff_matrix

        return CoordinateVector(diff_components + diff_base_vectors, self.system, self.point)

    @classmethod
    def from_expr(cls, expr: Expr) -> Self:
        combined = combine_coordinate_vectors(expr)

        if combined == 0:
            return combined

        if not isinstance(combined, cls):
            raise TypeError(f"Expected {cls.__name__} or zero, got {type(combined).__name__}")

        return combined


def combine_coordinate_vectors(expr: Expr) -> Expr:
    """
    Adds up coordinate vectors defined in the same coordinate system and at the same point within
    the given expression. Also takes into account the exact class of the vector, e.g.
    `QuantityCoordinateVector` is computed separately from pure `CoordinateVector`.

    Example:
    ========

    >>> from sympy import Symbol as SymSymbol
    >>> from symplyphysics.core.experimental.coordinate_systems import CylindricalCoordinateSystem, CoordinateVector
    >>> sys = CylindricalCoordinateSystem()
    >>> p = SymSymbol("P")
    >>> v1 = CoordinateVector([1, 2, 3], sys, p)
    >>> v2 = CoordinateVector([-1, -2, -3], sys, p)
    >>> assert combine_coordinate_vectors(v1 + v2) == 0
    """

    result = S.Zero

    mapping: dict[tuple[type[CoordinateVector], BaseCoordinateSystem, AppliedPoint | SymSymbol],
        ImmutableMatrix] = {}

    for term in into_terms(expr):
        vector, factor = split_factor(term)

        if not isinstance(vector, CoordinateVector):
            result += term
            continue

        mul = factor * vector.components

        args = type(vector), vector.system, vector.point

        total = mapping.get(args, ImmutableMatrix.zeros(3, 1))
        mapping[args] = total + mul

    for (cls, system, point), components in mapping.items():
        if components.is_zero_matrix:
            continue

        result += cls(components, system, point)

    return result


class QuantityCoordinateVector(CoordinateVector):
    _dimension: Dimension

    @property
    def dimension(self) -> Dimension:
        return self._dimension

    def __new__(
        cls,
        components: Iterable[Any],
        system: BaseCoordinateSystem,
        point: Optional[AppliedPoint | SymSymbol] = None,
    ) -> Expr:
        factors = []
        dimension = None

        for component in components:
            factor, dimension_ = collect_quantity_factor_and_dimension(component)
            factors.append(factor)

            if is_any_dimension(factor):
                continue

            if dimension is None:
                dimension = dimension_
            elif not dimsys_SI.equivalent_dims(dimension_, dimension):
                raise UnitsError(f"Expected {dimension}, got {dimension_}")

        dimension = dimension or dimensionless

        if dimsys_SI.is_dimensionless(dimension):
            components = factors
        else:
            components = [Quantity(factor, dimension=dimension) for factor in factors]

        obj = super().__new__(cls, components, system, point)

        if not isinstance(obj, CoordinateVector):
            return obj

        obj._dimension = dimension

        return obj

    # NOTE: Used in `symplyphysics.core.dimensions.collect_quantity` to avoid cyclic import
    def collect_quantity_factor_and_dimension(self) -> tuple[Expr, Dimension]:
        return self, self.dimension


__all__ = [
    "CoordinateVector",
    "QuantityCoordinateVector",
    "combine_coordinate_vectors",
]
