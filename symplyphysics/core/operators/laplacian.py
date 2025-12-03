from typing import Optional, Any

from sympy import Expr, S, diff

from ..miscellaneous import evaluate_or_global_fallback, sympify_expr
from ..vectors import is_vector_expr
from ..coordinate_systems import CoordinateScalar


class VectorLaplacian(Expr):  # pylint: disable=too-few-public-methods
    """
    The **Laplacian** is a differential operator given by the divergence of the gradient of a
    scalar function on Euclidean space.

    In an arbitrary curvilinear orthogonal coordinate system with base scalars `{q_i}` and Lamé
    coefficients `{h_i}`, the Laplacian of `f({q_i})` can be calculated as follows::

        Laplace(f) = (1 / J) * sum(diff((J / h_i^2) * diff(f, q_i), q_i), i)

    Here, `J = h_1 * h_2 * h_3` is the Jacobian determinant, and `diff` is the differentiation
    operator.
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

        system = scalar.system
        the_scalar = scalar.scalar
        base_scalars = system.base_scalars
        lame_coefficients = system.lame_coefficients(base_scalars)

        result = S.Zero

        for i in range(3):
            j, k = (i + 1) % 3, (i + 2) % 3

            ratio = lame_coefficients[j] * lame_coefficients[k] / lame_coefficients[i]

            q_i = base_scalars[i]

            result += diff(diff(the_scalar, q_i) * ratio, q_i)

        h1, h2, h3 = lame_coefficients
        result /= (h1 * h2 * h3)

        return CoordinateScalar(result, system, scalar.point)
