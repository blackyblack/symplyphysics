from __future__ import annotations
from typing import Any
from sympy import Expr, Basic, Idx, Sum


class IndexedSum(Expr):
    """
    Represents unevaluated Sum for expression with indexed variables.

    """

    def __new__(cls, function: Basic, *index_base: Idx, **assumptions: Any) -> IndexedSum:
        obj = Expr.__new__(cls, **assumptions)
        arglist = [function]
        if len(index_base) > 1:
            raise ValueError("Only one index is supported for Sum")
        arglist.extend(index_base)
        obj._args = tuple(arglist)
        return obj

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        return self.doit()._eval_nseries(x, n, logx, cdir)

    @property
    def index_base(self) -> Basic:
        return self._args[1]

    def doit(self, **hints: Any) -> Basic:
        # NOTE: This allows `IndexedSum` to be used in solvers, since `Sum` cannot be instantiated
        # without a bounded index.
        if len(self.index_base.args) < 2:
            return self

        inner_sum = Sum(self.args[0], self.index_base)
        return inner_sum.doit(**hints)
