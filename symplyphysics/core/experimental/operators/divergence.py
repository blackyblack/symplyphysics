from typing import Optional, Any

from sympy import Expr, S, diff, prod

from ..miscellaneous import evaluate_or_global_fallback
from ..vectors import is_vector_expr
from ..coordinate_systems import CoordinateScalar, CoordinateVector


class VectorDivergence(Expr):  # pylint: disable=too-few-public-methods
    """
    **Divergence** is a differential operator that operates on a vector field, producing a scalar
    field that gives the rate that the vector field alters the volume in an infinitesimal
    neighborhood at each point. In physical terms, it is the extent to which the vector field
    behaves like a source or a sink at each point.

    In an arbitrary curvilinear orthogonal coordinate system with base scalars `{q_i}` and Lamé
    coefficients `{h_i}`, the gradient of `F({q_i})` can be calculated as follows::

        div(F) = sum(diff(F_i * J / h_i, q_i), i) / J

    Here, `diff` is the differentiation operator, and `J = h_1 * h_2 * h_3` is the Jacobian
    determinant.

    **Links:**

    #. `Wikipedia <https://en.wikipedia.org/wiki/Divergence>`__.
    """

    def __new__(
        cls,
        vector: Any,
        *,
        evaluate: Optional[bool] = None,
    ) -> Expr:
        if vector == 0:
            return S.Zero

        if not is_vector_expr(vector):
            raise ValueError(f"Expected vector, got scalar: {vector}")

        evaluate = evaluate_or_global_fallback(evaluate)

        if not evaluate:
            obj = super().__new__(cls)  # pylint: disable=no-value-for-parameter
            obj._args = (vector,)

            return obj

        if not isinstance(vector, CoordinateVector):
            return S.Zero

        system = vector.system
        base_scalars = system.base_scalars
        lame_coefficients = system.lame_coefficients(base_scalars)

        result = S.Zero

        for i, (v_i, q_i) in enumerate(zip(vector.components, base_scalars)):
            h_jk = prod(lame_coefficients[j] for j in range(3) if j != i)

            result += diff(h_jk * v_i, q_i)

        h1, h2, h3 = lame_coefficients
        result /= (h1 * h2 * h3)

        return CoordinateScalar(result, system, vector.point)

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass
