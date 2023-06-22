from functools import reduce
import operator
from typing import Iterable
from sympy import Expr, Tuple


class SumArray(Expr):
    """
    Represents unevaluated Sum over array.

    """

    def __new__(cls, array):
        array = array if isinstance(array, (tuple, Tuple)) else (array,)
        obj = Expr.__new__(cls, *array)
        return obj

    def doit(self, **hints):
        return reduce(operator.add, self.args) if isinstance(self.args, Iterable) else self.args
