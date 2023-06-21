from functools import reduce
import operator
from typing import Iterable
from sympy import Expr


class SumArray(Expr):
    """
    Represents unevaluated Sum over array.

    """

    def __new__(cls, array):
        array = array if isinstance(array, tuple) else (array,)
        obj = Expr.__new__(cls, *array)
        return obj

    def doit(self, **hints):
        if len(self.args) == 0:
            return 0
        # Only one argument is expected in SumArray() and we always store
        # arguments as tuple
        return reduce(operator.add, self.args[0]) if isinstance(self.args[0], Iterable) else self.args[0]
