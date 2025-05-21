from typing import Optional, Any

from sympy import Expr, S, diff, prod

from ..miscellaneous import set_evaluate
from ..vectors import is_vector_expr
from ..coordinate_systems import CoordinateScalar, CoordinateVector


class VectorDivergence(Expr):  # pylint: disable=too-few-public-methods

    def __new__(
        cls,
        vector: Any,
        *,
        evaluate: Optional[bool] = None,
    ) -> Expr:
        if not is_vector_expr(vector):
            raise ValueError(f"Expected vector, got scalar: {vector}")

        if vector == 0 or not isinstance(vector, CoordinateVector):
            return S.Zero

        evaluate = set_evaluate(evaluate)

        if not evaluate:
            obj = super().__new__(cls)  # pylint: disable=no-value-for-parameter
            obj._args = (vector,)

            return obj

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
