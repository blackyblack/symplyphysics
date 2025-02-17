from typing import Any, Self
from sympy import Expr, Basic, Idx, Product


class IndexedProduct(Expr):
    """
    Represents unevaluated product for expression with indexed variables.

    """

    def __new__(cls, function: Basic, *index_base: Idx, **assumptions: Any) -> Self:
        obj = Expr.__new__(cls, **assumptions)
        arglist = [function]
        if len(index_base) > 1:
            raise ValueError("Only one index is supported for Product")
        arglist.extend(index_base)
        obj._args = tuple(arglist)
        return obj

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    @property
    def index_base(self) -> Basic:
        return self._args[1]

    def doit(self, **hints: Any) -> Basic:
        inner_sum = Product(self.args[0], self.index_base)
        return inner_sum.doit(**hints)
