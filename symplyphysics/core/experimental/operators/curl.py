from typing import Optional, Any
from itertools import permutations

from sympy import Expr, S, diff, Matrix
from sympy.combinatorics.permutations import Permutation

from ..miscellaneous import set_evaluate
from ..vectors import is_vector_expr, VectorExpr
from ..coordinate_systems import CoordinateVector


class VectorCurl(VectorExpr):  # pylint: disable=too-few-public-methods

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

        results: Matrix = Matrix.zeros(3, 1)

        for i, j, k in permutations(range(3)):
            e_ijk = Permutation([i, j, k]).signature()

            h_i = lame_coefficients[i]
            h_k = lame_coefficients[k]
            v_k = vector.components[k]
            q_j = base_scalars[j]

            results[i] += e_ijk * h_i * diff(h_k * v_k, q_j)

        h1, h2, h3 = lame_coefficients
        results /= h1 * h2 * h3

        return CoordinateVector(results, system, vector.point)
