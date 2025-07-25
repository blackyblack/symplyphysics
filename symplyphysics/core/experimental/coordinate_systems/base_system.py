from __future__ import annotations

from typing import Optional, Sequence, Iterable

from sympy import Expr, sqrt, Matrix, true, simplify, Symbol as SymSymbol, ImmutableMatrix, Basic
from sympy.logic.boolalg import Boolean

from symplyphysics.core.symbols.symbols import Symbol, clone_as_function, Function

from ..miscellaneous import const

_base_coordinate_system_cache: dict[tuple[type[BaseCoordinateSystem], Optional[tuple[Symbol, Symbol,
    Symbol]], Optional[tuple[Function, function, function]]], BaseCoordinateSystem] = {}


class BaseCoordinateSystem(Basic):
    _base_scalars: tuple[Symbol, Symbol, Symbol]

    _cartesian_derivative_matrix: ImmutableMatrix
    """See `BaseCoordinateSystem.cartesian_derivative_matrix`."""

    _lame_coefficients: tuple[Expr, Expr, Expr]
    """See `BaseCoordinateSystem.lame_coefficients`."""

    _base_vector_matrix: ImmutableMatrix
    """See `BaseCoordinateSystem.base_vector_matrix`."""

    _wrt: Symbol
    _base_scalar_functions: tuple[Function, Function, Function]
    _diff_base_vector_matrix: ImmutableMatrix
    """See `BaseCoordinateSystem.diff_base_vector_matrix`."""

    @property
    def base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return self._base_scalars

    @property
    def base_scalar_functions(self) -> tuple[Function, Function, Function]:
        return self._base_scalar_functions

    def __repr__(self) -> str:
        return type(self).__qualname__

    def __new__(
        cls,
        base_scalars: Optional[tuple[Symbol, Symbol, Symbol]] = None,
        base_scalar_functions: Optional[tuple[Function, Function, Function]] = None,
    ) -> BaseCoordinateSystem:
        if base_scalars:
            base_scalars = tuple(base_scalars)

        if base_scalar_functions:
            base_scalar_functions = tuple(base_scalar_functions)

        cached = _base_coordinate_system_cache.get((cls, base_scalars, base_scalar_functions))

        if cached is not None:
            return cached

        obj = super().__new__(cls)  # pylint: disable=no-value-for-parameter
        _base_coordinate_system_cache[cls, base_scalars, base_scalar_functions] = obj

        if not base_scalars:
            base_scalars = obj.generate_base_scalars()

        a, b, c = base_scalars
        obj._base_scalars = a, b, c

        x, y, z = obj.cartesian_transform(base_scalars)

        obj._cartesian_derivative_matrix = ImmutableMatrix([
            [x.diff(a), y.diff(a), z.diff(a)],
            [x.diff(b), y.diff(b), z.diff(b)],
            [x.diff(c), y.diff(c), z.diff(c)],
        ])

        obj._lame_coefficients = obj.generate_lame_coefficients()

        obj._base_vector_matrix = ImmutableMatrix([
            simplify(obj._cartesian_derivative_matrix.row(i_row) / obj._lame_coefficients[i_row])
            for i_row in range(3)
        ])

        obj._wrt = Symbol("t", real=True)

        if not base_scalar_functions:
            base_scalar_functions = (
                clone_as_function(a, [obj._wrt]),
                clone_as_function(b, [obj._wrt]),
                clone_as_function(c, [obj._wrt]),
            )

        obj._base_scalar_functions = base_scalar_functions

        applied_base_scalars = [f(obj._wrt) for f in obj._base_scalar_functions]

        diff_matrix = simplify(obj.base_vector_matrix(applied_base_scalars).diff(obj._wrt))
        inv_matrix = obj.inverse_base_vector_matrix(applied_base_scalars)

        obj._diff_base_vector_matrix = ImmutableMatrix(simplify(diff_matrix * inv_matrix))

        return obj

    def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        raise NotImplementedError

    def cartesian_transform(
        self,
        base_scalars: Optional[Sequence[Expr]] = None,
    ) -> Sequence[Expr]:
        """
        Defines how Cartesian scalars `x, y, z` can be expressed using base scalars `{q_i}` of the
        system in question.
        """

        raise NotImplementedError

    def inverse_cartesian_transform(self, cartesian_scalars: Sequence[Expr]) -> Sequence[Expr]:
        """
        Defined how base scalars `{q_i}` of the system in question can be expressed using Cartesian
        scalars `x, y, z`.
        """

        raise NotImplementedError

    def assumption(self, base_scalars: Optional[Sequence[Expr]] = None) -> Boolean:  # pylint: disable=unused-argument
        """
        Assumption(s) about base scalars in the coordinate system. For example, in spherical
        coordinates `r ≥ 0` and `sin(θ) ≥ 0` where `r` is the radius and `θ` is the polar angle.
        """

        return true

    def base_scalar_subs(self, base_scalars: Optional[Sequence[Expr]] = None) -> dict[Symbol, Expr]:
        if not base_scalars:
            return {}

        return {old: new for old, new in zip(self.base_scalars, base_scalars) if old != new}

    def cartesian_derivative_matrix(self, base_scalars: Optional[Sequence[Expr]] = None) -> Matrix:
        """
        Matrix whose element at row `i` and column `j` is `d(x_j)/d(q_i)`, where `{q_i} = {q_0,
        q_1, q_2}` are the base scalars of the system in question, and `{x_j} = {x, y, z}` are
        Cartesian scalars being functions of `{q_i}`.
        """
        base_scalars = base_scalars or self.base_scalars

        matrix = Matrix(self._cartesian_derivative_matrix)

        subs = self.base_scalar_subs(base_scalars)
        matrix = matrix.subs(subs)

        return simplify(matrix)

    def generate_lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        assumption = self.assumption(self.base_scalars)

        lame_coefficients = []

        for i_row in range(3):
            row = self._cartesian_derivative_matrix.row(i_row)
            lame_coefficient = sqrt(sum(elem**2 for elem in row)).refine(assumption).simplify()
            lame_coefficients.append(lame_coefficient)

        h1, h2, h3 = lame_coefficients

        return h1, h2, h3

    def lame_coefficients(self, base_scalars: Optional[Sequence[Expr]] = None) -> Sequence[Expr]:
        """
        An `i`-th Lamé coefficient is the length of the derivative of the position vector `r` with
        respect to the `i`-th base scalar of the system in question.
        """

        subs = self.base_scalar_subs(base_scalars)

        return tuple(
            lame_coefficient.subs(subs).simplify() for lame_coefficient in self._lame_coefficients)

    def base_vector_matrix(
        self,
        base_scalars: Optional[Sequence[Expr]] = None,
    ) -> Matrix:
        """
        Matrix whose rows represent the components of the base vectors of the system in question in
        terms of Cartesian base vectors, i.e. `A[i, j] = dot(Q_i, X_j)` where `A` is the matrix,
        `Q_i` is the `i`-th base vector of the given system, `X_j` is the `j`-th Cartesian base
        vector, and `dot` is the dot product.
        """

        matrix = Matrix(self._base_vector_matrix)

        subs = self.base_scalar_subs(base_scalars)
        matrix = matrix.subs(subs)

        return simplify(matrix)

    def inverse_base_vector_matrix(
        self,
        base_scalars: Optional[Sequence[Expr]] = None,
    ) -> Matrix:
        """
        Matrix whose rows represent the components of the Cartesian base vectors in terms of the
        base vectors of the system in question. It is exactly the transpose of the
        `base_vector_matrix`.
        """

        return self.base_vector_matrix(base_scalars).T

    def diff_base_vector_matrix(
        self,
        wrt: SymSymbol,
        base_scalars: Optional[Sequence[Expr]] = None,
    ) -> Matrix:
        """
        Returns a matrix whose element at row `i` and column `j` represents the `j`-th component
        of the derivative of the `i`-th base vector of the given system with respect to `wrt`, i.e.
        `A[i, j] = dot(diff(Q_i, wrt), Q_j)` where `A` is the matrix, `Q_i` is the `i`-th base
        vector of the given system, and `dot` is the dot product.
        """
        base_scalars = base_scalars or self.base_scalars

        diff_matrix = self._diff_base_vector_matrix.subs(self._wrt, wrt)

        matrix = Matrix.zeros(3)

        for i_row in range(diff_matrix.rows):
            for i_col in range(diff_matrix.cols):
                elem = diff_matrix[i_row, i_col]

                for func_base_scalar, base_scalar in zip(self._base_scalar_functions, base_scalars):
                    elem = elem.replace(func_base_scalar, const(base_scalar))

                matrix[i_row, i_col] = elem

        return simplify(matrix)

    def extract_position_vector_components(self, coordinates: Iterable[Expr]) -> Sequence[Expr]:
        """
        Convert the coordinates of a point into the components of the position vector that starts
        at the origin of the coordinate system and ends in that point.
        """

        raise NotImplementedError


__all__ = ["BaseCoordinateSystem"]
