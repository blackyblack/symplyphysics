from __future__ import annotations
from typing import Any
from sympy import Expr, Symbol as SymSymbol
from ..dimensions import Dimension, collect_dimension


class Symbolic(SymSymbol):  # type: ignore[misc]  # pylint: disable=too-many-ancestors
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

    def __new__(
        cls,
        expr: Expr,
        *,
        wrap_code: bool = False,
        wrap_latex: bool = False,
        **assumptions: Any,
    ) -> Symbolic:
        cls_name = cls.__name__
        inner = str(expr)
        display_name = f"{cls_name}({inner})"

        obj = super().__new__(cls, display_name, **assumptions)
        return obj  # type: ignore[no-any-return]

    def __init__(
        self,
        expr: Expr,
        *,
        wrap_code: bool = False,
        wrap_latex: bool = False,
    ) -> None:
        self.factor = expr
        self.dimension = collect_dimension(expr)
        self.wrap_code = wrap_code
        self.wrap_latex = wrap_latex


class Average(Symbolic):  # pylint: disable=too-many-ancestors
    """
    Represents the average value of the argument :math:`x`.

    Code:
        :code:`avg(x)`

    Latex:
        :math:`\\langle x \\rangle`

    **Notes:**

    #. The `wrap_code` and `wrap_latex` fields are unused.
    """


class FiniteDifference(Symbolic):  # pylint: disable=too-many-ancestors
    """
    Represents the finite difference of the argument :math:`x`, e.g.
    :math:`\\Delta x = x_2 - x_1`.

    Code:
        :code:`Delta(x)`

    Latex:
        :math:`\\Delta x`

    **Notes:**

    #. The `wrap_code` field is unused.
    """


class ExactDifferential(Symbolic):  # pylint: disable=too-many-ancestors
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


class InexactDifferential(Symbolic):  # pylint: disable=too-many-ancestors
    """
    A differential is said to be **inexact** (or **imperfect**) when its integral is
    path dependent.

    Code:
        :code:`delta(Gamma)`

    Latex:
        :math:`\\delta \\Gamma`

    **Links:**

    #. `Wikipedia <https://en.wikipedia.org/wiki/Inexact_differential#>`__.

    **Notes:**

    #. The `wrap_code` field is unused.
    """
