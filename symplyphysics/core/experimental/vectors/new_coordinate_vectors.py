from __future__ import annotations

from typing import Sequence, Optional, Iterable, SupportsFloat, Sized, Any
from sympy import (ImmutableMatrix, Expr, sqrt, Symbol as SymSymbol, Basic, Derivative as
    SymDerivative, Atom)
from sympy.printing.printer import Printer
from sympy.matrices.dense import DenseMatrix

from symplyphysics.core.expr_comparisons import expr_equals

from ..miscellaneous import sympify_expr
from ..coordinate_systems.new_coordinate_systems import BaseCoordinateSystem, CartesianCoordinateSystem
from . import VectorExpr, is_vector_expr, into_terms, split_factor, VectorCross


def _prepare(coordinates: Iterable[SupportsFloat], system: BaseCoordinateSystem) -> dict[SymSymbol,
    Expr]:
    return {
        scalar: sympify_expr(coordinate)
        for scalar, coordinate in zip(system.base_scalars, coordinates)
    }


class AppliedPoint(Atom):
    """
    An `AppliedPoint` corresponds to a point in (3D) space whose coordinates are defined within a
    certain coordinate system.
    """

    _coordinates: dict[SymSymbol, Expr]
    _system: BaseCoordinateSystem

    _iterable = False

    @property
    def coordinates(self) -> dict[SymSymbol, Expr]:
        return self._coordinates

    @property
    def system(self) -> BaseCoordinateSystem:
        return self._system

    def __getitem__(self, base_scalar: SymSymbol) -> Expr:
        return self.coordinates[base_scalar]

    def __new__(
        cls,
        coordinates: Iterable[Any],
        system: BaseCoordinateSystem,
    ) -> AppliedPoint:
        return super().__new__(cls)  # pylint: disable=no-value-for-parameter

    def __init__(
        self,
        coordinates: Iterable[Any],
        system: BaseCoordinateSystem,
    ):
        super().__init__()

        if isinstance(coordinates, Sized):
            n = len(coordinates)  # can't extract out of if-block, mypy complains otherwise
        else:
            coordinates = tuple(coordinates)
            n = len(coordinates)

        if n != 3:
            raise ValueError(f"The point must have all 3 coordinates defined, got {n}.")

        self._coordinates = _prepare(coordinates, system)
        self._system = system

    def _sympystr(self, p: Printer) -> str:
        system_name = type(self.system).__name__.removesuffix("CoordinateSystem")
        point_name = f"{system_name}Point"

        coordinates = ", ".join(
            f"{p.doprint(s)} = {p.doprint(c)}" for s, c in self.coordinates.items())

        return f"{point_name}({coordinates})"

    def _hashable_content(self) -> tuple[Any, ...]:
        return (tuple(self.coordinates.items()),)

    def equals(self, other: AppliedPoint) -> bool:
        for base_scalar, coordinate in self.coordinates.items():
            if base_scalar not in other.coordinates:
                return False

            if not expr_equals(coordinate, other[base_scalar]):
                return False

        return True


# Used as a common point for all Cartesian vectors
_CartesianVectorPoint = SymSymbol("P")


class CoordinateVector(VectorExpr):
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
    def free_symbols(self) -> set[Expr]:
        return self.components.free_symbols  # type: ignore[no-any-return]

    def _hashable_content(self) -> tuple[Basic, ...]:
        return self.components, self.system, self.point

    def __init__(
        self,
        components: Sequence[Expr],
        system: BaseCoordinateSystem,
        point: Optional[AppliedPoint | SymSymbol] = None,
    ) -> None:
        super().__init__()

        if isinstance(components, DenseMatrix) and 1 not in components.shape:
            rows, cols = components.shape
            message = f"Expected a row or column vector, got a {rows} by {cols} matrix"
            raise ValueError(message)

        if isinstance(system, CartesianCoordinateSystem):
            # Cartesian vectors are independent of their point of application.
            point = _CartesianVectorPoint
        elif point is None:
            message = "The point of application must be defined for vectors in non-Cartesian systems"
            raise ValueError(message)

        if isinstance(point, AppliedPoint) and system != point.system:
            message = "The system of the vector and the point of its application must coincide"
            raise ValueError(message)

        components = [sympify_expr(component) for component in components]
        self._components = ImmutableMatrix(components)

        self._system = system
        self._point = point

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


def combine_coordinate_vectors(expr: Expr) -> Expr:
    non_coordinate = []
    system_to_components: dict[tuple[BaseCoordinateSystem, Optional[AppliedPoint | SymSymbol]],
        ImmutableMatrix] = {}

    for cross in expr.atoms(VectorCross):
        lhs, rhs = cross.args
        lhs = combine_coordinate_vectors(lhs)
        rhs = combine_coordinate_vectors(rhs)
        expr = expr.subs(cross, VectorCross(lhs, rhs))

    for term in into_terms(expr):
        vector, factor = split_factor(term)

        if not isinstance(vector, CoordinateVector):
            non_coordinate.append(term)
            continue

        total = system_to_components.get((vector.system, vector.point), ImmutableMatrix.zeros(3, 1))
        mul = factor * vector.components
        system_to_components[vector.system, vector.point] = total + mul

    result = sum(non_coordinate)

    for (system, point), components in system_to_components.items():
        if components.is_zero_matrix:
            continue

        result += CoordinateVector(components, system, point)

    return result
