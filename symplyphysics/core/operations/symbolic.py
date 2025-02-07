from __future__ import annotations
from typing import Any
from sympy import Expr
from symplyphysics.core.dimensions import Dimension, collect_factor_and_dimension


class Symbolic(Expr):
    """
    This class is intended to be subclassed for custom code printing.
    """

    factor: Expr
    """
    Argument of the wrapper.
    """

    dimension: Dimension
    """
    Dimension of the argument.
    """

    wrap_code: bool
    """
    Denotes whether the code printer should wrap the argument in parentheses.
    """

    wrap_latex: bool
    """
    Denotes whether the latex printer should wrap the argument in parentheses.
    """

    def __new__(cls, *args: Any, **kwargs: Any) -> Symbolic:
        return super().__new__(cls)

    def __init__(self, expr: Expr, wrap_code: bool = False, wrap_latex: bool = False) -> None:
        factor, dimension = collect_factor_and_dimension(expr)
        self.factor = factor
        self.dimension = dimension
        self.wrap_code = wrap_code
        self.wrap_latex = wrap_latex

    def _eval_nseries(self, _x: Any, _n: Any, _logx: Any, _cdir: Any) -> Any:
        pass


class Average(Symbolic):
    """
    Represents the average value of the argument :math:`x`.

    Code:
        :code:`avg(x)`

    Latex:
        :math:`\\langle x \\rangle`
    """


class FiniteDifference(Symbolic):
    """
    Represents the finite difference of the argument :math:`x`, e.g.
    :math:`\\Delta x = x_2 - x_1`.

    Code:
        :code:`Delta(x)`

    Latex:
        :math:`\\Delta x`
    """


class ExactDifferential(Symbolic):
    """
    A differential is said to be **exact** (or **perfect**) if it is equal to the
    general differential :math:`d \\Gamma` for some differentiable function
    :math:`\\Gamma` in an orthogonal coordinate system. Equivalently, its integral is
    path independent.

    Code:
        :code:`d(Gamma)`

    Latex:
        :math:`d \\Gamma`

    **Links:**
    
    #. `Wikipedia <https://en.wikipedia.org/wiki/Exact_differential#>`__.
    """


class InexactDifferential(Symbolic):
    """
    A differential is said to be **inexact** (or **imperfect**) when its integral is
    path dependent.

    Code:
        :code:`delta(Gamma)`

    Latex:
        :math:`\\delta \\Gamma`

    **Links:**

    #. `Wikipedia <https://en.wikipedia.org/wiki/Inexact_differential#>`__.
    """
