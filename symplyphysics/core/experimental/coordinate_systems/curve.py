from __future__ import annotations

from sympy import Basic, Symbol as SymSymbol, Expr, Integral

from .point import AppliedPoint


class Curve(Basic):
    """
    Topologically, a **curve** is specified by a continuous function :math:`\\gamma: I \\to C`
    from an interval :math:`I` of real numbers into a topological space :math:`X`.
    """

    @property
    def parameter(self) -> SymSymbol:
        """A symbol that parametrized the curve."""

        return self.args[0]

    @property
    def parametrization(self) -> AppliedPoint:
        """Representation of the curve as a function of the curve parameter."""

        return self.args[1]

    def __new__(
        cls,
        parameter: SymSymbol,
        parametrization: AppliedPoint,
    ) -> Curve:
        return super().__new__(cls, parameter, parametrization)  # pylint: disable=too-many-function-args

    def to_integral(self, expr: Expr, bounds: tuple[Expr, Expr]) -> Integral:
        lo, hi = bounds
        return Integral(expr, (self.parameter, lo, hi))
