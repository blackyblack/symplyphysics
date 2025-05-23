from typing import Optional, Any

from sympy import Expr, S, ImmutableMatrix, diff

from ..miscellaneous import evaluate_or_global_fallback, sympify_expr
from ..vectors import VectorExpr, is_vector_expr
from ..coordinate_systems import CoordinateScalar, CoordinateVector


class VectorGradient(VectorExpr):  # pylint: disable=too-few-public-methods
    """
    Assuming `f` is a scalar-values differentiable function of several variables, its **gradient**
    is a vector field (i.e. a vector-valued function) `grad(f)` whose value at a point `P` gives
    the direction and the rate of fastest increase.

    In an arbitrary curvilinear orthogonal coordinate system with base scalars `{q_i}` and Lamé
    coefficients `{h_i}`, the gradient of `f({q_i})` can be calculated as follows::

        grad(f) = sum(diff(f, q_i) / h_i * e_i, i)

    Here, `diff` is the differentiation operator, and `e_i` is the `i`-th unit base vector.

    **Links:**

    #. `Wikipedia <https://en.wikipedia.org/wiki/Gradient>`__.
    """

    def __new__(
        cls,
        scalar: Any,
        *,
        evaluate: Optional[bool] = None,
    ) -> Expr:
        if scalar == 0:
            return S.Zero

        scalar = sympify_expr(scalar)

        if is_vector_expr(scalar):
            raise ValueError(f"Expected scalar, got vector: {scalar}")

        evaluate = evaluate_or_global_fallback(evaluate)

        if not evaluate:
            obj = super().__new__(cls)  # pylint: disable=no-value-for-parameter
            obj._args = (scalar,)

            return obj

        if not isinstance(scalar, CoordinateScalar):
            return S.Zero

        the_scalar = scalar.scalar
        system = scalar.system
        base_scalars = system.base_scalars
        lame_coefficients = system.lame_coefficients(base_scalars)

        components = ImmutableMatrix(
            [diff(the_scalar, q) / h for q, h in zip(base_scalars, lame_coefficients)])

        return CoordinateVector(components, system, scalar.point)


__all__ = ["VectorGradient"]
