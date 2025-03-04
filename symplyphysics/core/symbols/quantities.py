from __future__ import annotations

from functools import partial
from typing import Any, Optional, Sequence, Any
from sympy import S, Expr, sympify, Abs
from sympy.physics.units import Dimension, Quantity as SymQuantity
from sympy.physics.units.systems.si import SI
from sympy.multipledispatch import dispatch
from sympy.printing.printer import Printer

from .symbols import DimensionSymbol, next_name
from ..dimensions.collect_quantity import collect_quantity_factor_and_dimension
from ..dimensions import dimension_to_si_unit


class Quantity(DimensionSymbol, SymQuantity):  # type: ignore[misc]  # pylint: disable=too-many-ancestors

    # pylint: disable-next=signature-differs
    def __new__(cls,
        _expr: Any = S.One,
        *,
        display_symbol: Optional[str] = None,
        display_latex: Optional[str] = None,
        **assumptions: Any) -> Quantity:
        name = next_name("QTY")
        # Latex symbol is set in SymPy Quantity, not in DimensionSymbol, due to Latex printer
        # specifics
        display_symbol = display_symbol or name
        display_latex = display_latex or display_symbol
        obj = SymQuantity.__new__(cls, name, None, display_latex, None, None, None, False,
            **assumptions)
        return obj  # type: ignore[no-any-return]

    def __init__(self,
        expr: Any = S.One,
        *,
        display_symbol: Optional[str] = None,
        display_latex: Optional[str] = None,
        dimension: Optional[Dimension] = None) -> None:
        (scale, dimension_) = collect_quantity_factor_and_dimension(expr)
        try:
            _ = complex(scale)  # if this fails, then ``scale`` contains a symbolic sub-expression
        except Exception as e:
            raise ValueError(
                f"Argument '{expr}' to function 'Quantity()' should "
                f"be an expression made of numbers and quantities.",) from e

        dimension = dimension or dimension_
        display_symbol = display_symbol or str(self.name)
        super().__init__(display_symbol, dimension, display_latex=display_latex)
        SI.set_quantity_dimension(self, dimension)
        SI.set_quantity_scale_factor(self, scale)

    # This is required for integration to work properly
    @property
    def func(self) -> partial[Quantity]:
        return partial(Quantity.identity, self)

    def identity(self, *_args: Any) -> Quantity:
        return self

    def _eval_is_positive(self) -> bool:
        # NOTE: returns False for complex values, see https://github.com/blackyblack/symplyphysics/blob/3e7e05b9837c70bb23d36202b9e958b739cd36bc/test/electricity/circuits/transmission_lines/transmission_matrix_lossy_transmission_line_test.py#L23
        try:
            return scale_factor(self) >= 0
        except TypeError:
            return False

    def _eval_Abs(self) -> Quantity:
        return self.__class__(Abs(self.scale_factor), dimension=self.dimension)

    def _sympystr(self, p: Printer) -> str:
        if "QTY" not in self.display_name:
            return self.display_name

        si_unit = dimension_to_si_unit(self.dimension)

        si_value = self.convert_to(si_unit) / si_unit

        qty: SymQuantity

        for qty in si_value.atoms(SymQuantity):
            si_value = si_value.subs(qty, 1)
        si_value = si_value.n(3)

        for qty in si_unit.atoms(SymQuantity):
            abbrev = qty.abbrev
            if abbrev:
                si_unit = si_unit.subs(qty, abbrev)

        return str(p.doprint(si_value * si_unit))


# Allows for some SymPy comparisons, eg Piecewise function
@dispatch(Quantity, Quantity)  # type: ignore[misc]
def _eval_is_ge(lhs: Quantity, rhs: Quantity) -> bool:
    return scale_factor(lhs) >= scale_factor(rhs)


def subs_list(input_: Sequence[Expr | float], subs_: dict[Expr, Quantity]) -> Sequence[Quantity]:
    return [Quantity(sympify(c).subs(subs_)) for c in input_]


def scale_factor(quantity_: Quantity | float) -> float:
    if isinstance(quantity_, Quantity):
        return float(quantity_.scale_factor)

    return quantity_
