from __future__ import annotations
from typing import Any
from sympy import Expr, Basic, Idx, Product


class IndexedProduct(Expr):  # type: ignore[misc]
    """
    Represents unevaluated product for expression with indexed variables.

    """

    def __new__(cls, function: Basic, *index_base: Idx, **assumptions: Any) -> IndexedProduct:
        obj = Expr.__new__(cls, **assumptions)
        arglist = [function]
        if len(index_base) > 1:
            raise ValueError("Only one index is supported for Product")
        arglist.extend(index_base)
        obj._args = tuple(arglist)
        return obj  # type: ignore[no-any-return]

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    @property
    def index_base(self) -> Basic:
        return self._args[1]

    def doit(self, **hints: Any) -> Basic:
        inner_sum = Product(self.args[0], self.index_base)
        return inner_sum.doit(**hints)
