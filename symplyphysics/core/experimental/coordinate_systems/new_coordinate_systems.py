from __future__ import annotations

from typing import Optional, Sequence
from sympy import (Expr, sqrt, Matrix, true, sin, cos, tan, Q, simplify, atan2, Symbol as SymSymbol,
    S, ImmutableMatrix, Basic)
from sympy.logic.boolalg import Boolean
from symplyphysics import Symbol, units, clone_as_function, Function

from ..miscellaneous import cacheit


class BaseCoordinateSystem(Basic):  # type: ignore[misc]
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

    def _hashable_content(self) -> tuple[Basic, ...]:
        return ()

    def __repr__(self) -> str:
        return type(self).__qualname__

    @cacheit
    def __new__(
        cls,
        base_scalars: Optional[tuple[Symbol, Symbol, Symbol]] = None,
        base_scalar_functions: Optional[tuple[Function, Function, Function]] = None,
    ) -> BaseCoordinateSystem:
        obj = super().__new__(cls)

        if not base_scalars:
            base_scalars = obj.generate_base_scalars()

        q1, q2, q3 = base_scalars
        obj._base_scalars = q1, q2, q3

        a, b, c = base_scalars
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

        return obj  # type: ignore[no-any-return]

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
                    elem = elem.replace(func_base_scalar, lambda _: base_scalar)  # pylint: disable=cell-var-from-loop

                matrix[i_row, i_col] = elem

        return simplify(matrix)


class CartesianCoordinateSystem(BaseCoordinateSystem):

    def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return (
            Symbol("x", units.length, real=True),
            Symbol("y", units.length, real=True),
            Symbol("z", units.length, real=True),
        )

    def cartesian_transform(self, base_scalars: Optional[Sequence[Expr]] = None) -> Sequence[Expr]:
        return base_scalars or self.base_scalars

    def inverse_cartesian_transform(self, cartesian_scalars: Sequence[Expr]) -> Sequence[Expr]:
        return cartesian_scalars

    def generate_lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        return S.One, S.One, S.One

    def cartesian_derivative_matrix(
        self,
        _base_scalars: Optional[Sequence[Expr]] = None,
    ) -> Matrix:
        return Matrix.eye(3)

    def base_vector_matrix(
        self,
        _base_scalars: Optional[Sequence[Expr]] = None,
    ) -> Matrix:
        return Matrix.eye(3)

    def diff_base_vector_matrix(
        self,
        _wrt: SymSymbol,
        _base_scalars: Optional[Sequence[Expr]] = None,
    ) -> Matrix:
        return Matrix.zeros(3)


class CylindricalCoordinateSystem(BaseCoordinateSystem):

    def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return (
            Symbol("rho", units.length, display_latex="\\rho", nonnegative=True),
            Symbol("phi", display_latex="\\varphi", real=True),
            Symbol("z", units.length, real=True),
        )

    def cartesian_transform(self, base_scalars: Optional[Sequence[Expr]] = None) -> Sequence[Expr]:
        rho, phi, z = base_scalars or self.base_scalars

        x = rho * cos(phi)
        y = rho * sin(phi)

        return x, y, simplify(z)

    def inverse_cartesian_transform(self, cartesian_scalars: Sequence[Expr]) -> Sequence[Expr]:
        x, y, z = cartesian_scalars

        rho = sqrt(x**2 + y**2)
        phi = atan2(y, x)

        return rho, phi, z

    def assumption(self, base_scalars: Optional[Sequence[Expr]] = None) -> Boolean:
        rho = (base_scalars or self.base_scalars)[0]

        return Q.positive(rho)  # pylint: disable=too-many-function-args

    def generate_lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        rho, _, _ = self.base_scalars

        h_rho = S.One
        h_phi = rho
        h_z = S.One

        return h_rho, h_phi, h_z


class SphericalCoordinateSystem(BaseCoordinateSystem):

    def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return (
            Symbol("r", units.length, nonnegative=True),
            Symbol("theta", display_latex="\\theta", nonnegative=True),
            Symbol("phi", display_latex="\\varphi", real=True),
        )

    def cartesian_transform(self, base_scalars: Optional[Sequence[Expr]] = None) -> Sequence[Expr]:
        r, theta, phi = base_scalars or self.base_scalars

        x = r * sin(theta) * cos(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(theta)

        return x, y, z

    def inverse_cartesian_transform(self, cartesian_scalars: Sequence[Expr]) -> Sequence[Expr]:
        x, y, z = cartesian_scalars

        r = sqrt(x**2 + y**2 + z**2)
        theta = atan2(sqrt(x**2 + y**2), z)
        phi = atan2(y, x)

        return r, theta, phi

    def assumption(self, base_scalars: Optional[Sequence[Expr]] = None) -> Boolean:
        r, theta, *_ = base_scalars or self.base_scalars

        # pylint: disable-next=too-many-function-args
        return Q.positive(r) & Q.positive(sin(theta)) & Q.positive(tan(theta))

    def generate_lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        r, theta, _ = self.base_scalars

        h_r = S.One
        h_theta = r
        h_phi = r * sin(theta)

        return h_r, h_theta, h_phi
