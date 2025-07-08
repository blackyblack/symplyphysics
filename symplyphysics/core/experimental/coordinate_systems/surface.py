from __future__ import annotations

from typing import Iterable

from sympy import Basic, Symbol as SymSymbol, Expr, Integral

from .point import AppliedPoint


class Surface(Basic):
    """
    Topologically, a **surface** is TODO
    """

    @property
    def parameters(self) -> tuple[SymSymbol, SymSymbol]:
        """The symbols that parametrized the surface."""

        return self.args[0]

    @property
    def parametrization(self) -> AppliedPoint:
        """Representation of the surface as a function of the parameters."""

        return self.args[1]

    def __new__(
        cls,
        parameters: Iterable[SymSymbol] | tuple[SymSymbol, SymSymbol],
        parametrization: AppliedPoint,
    ) -> Surface:
        p1, p2 = parameters

        return super().__new__(cls, (p1, p2), parametrization)  # pylint: disable=too-many-function-args

    def to_integral(
        self,
        expr: Expr,
        first_bounds: tuple[Expr, Expr],
        second_bounds: tuple[Expr, Expr],
    ) -> Integral:
        p1, p2 = self.parameters

        lo1, hi1 = first_bounds
        lo2, hi2 = second_bounds

        return Integral(expr, (p1, lo1, hi1), (p2, lo2, hi2))
