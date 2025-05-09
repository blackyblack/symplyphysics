from typing import Sequence, Optional
from sympy import ImmutableMatrix, Expr, sqrt, Symbol as SymSymbol, Basic, Derivative as SymDerivative
from sympy.matrices.dense import DenseMatrix

from ..miscellaneous import sympify_expr
from ..coordinate_systems.new_coordinate_systems import BaseCoordinateSystem, CartesianCoordinateSystem
from . import VectorExpr, is_vector_expr, into_terms, split_factor, VectorCross


class CoordinateVector(VectorExpr):
    _components: ImmutableMatrix
    _system: BaseCoordinateSystem

    is_number = False
    is_atomic_vector = True

    @property
    def components(self) -> ImmutableMatrix:
        return self._components

    @property
    def system(self) -> BaseCoordinateSystem:
        return self._system

    @property
    def free_symbols(self) -> set[Expr]:
        return self.components.free_symbols  # type: ignore[no-any-return]

    def _hashable_content(self) -> tuple[Basic, ...]:
        return self.components, self.system

    def __init__(self, components: Sequence[Expr], system: BaseCoordinateSystem) -> None:
        if isinstance(components, DenseMatrix):
            if 1 not in components.shape:
                rows, cols = components.shape
                message = f"Expected a row or column vector, got a {rows} by {cols} matrix"
                raise ValueError(message)

        components = [sympify_expr(component) for component in components]
        self._components = ImmutableMatrix(components)

        self._system = system

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

        if type(lhs.system) is not type(rhs.system):
            return None

        components = lhs.components.cross(rhs.components)

        return CoordinateVector(components, lhs.system)

    def _eval_derivative(self, symbol: SymSymbol) -> VectorExpr:
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

        return CoordinateVector(diff_components + diff_base_vectors, self.system)


def combine_coordinate_vectors(expr: Expr) -> Expr:
    non_coordinate = []
    system_to_components: dict[BaseCoordinateSystem, ImmutableMatrix] = {}

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

        total = system_to_components.get(vector.system, ImmutableMatrix.zeros(3, 1))
        mul = factor * vector.components
        system_to_components[vector.system] = total + mul

    result = sum(non_coordinate)

    for system, components in system_to_components.items():
        if components.is_zero_matrix:
            continue

        result += CoordinateVector(components, system)

    return result
