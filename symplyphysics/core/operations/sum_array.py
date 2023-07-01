from __future__ import annotations
from functools import reduce
import operator
from typing import Any
from sympy import Expr, Tuple


class SumArray(Expr):
    """
    Represents unevaluated Sum over array.

    """

    def __new__(cls, array: Expr | Tuple | tuple[Expr, ...]) -> SumArray:
        array_unpacked = array if isinstance(array, (tuple, Tuple)) else Tuple(array)
        obj = Expr.__new__(cls, *array_unpacked)
        return obj

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def doit(self, **_hints: Any) -> Expr:
        return reduce(operator.add, self.args)
