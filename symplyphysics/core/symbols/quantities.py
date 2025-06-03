from __future__ import annotations

from functools import partial
from typing import Any, Optional, Sequence, SupportsFloat
from sympy import S, Expr, sympify, Abs
from sympy.physics.units import Dimension, Quantity as SymQuantity
from sympy.physics.units.systems.si import SI
from sympy.multipledispatch import dispatch
from sympy.printing.printer import Printer

from .symbols import DimensionSymbol, next_name
from ..dimensions.collect_quantity import collect_quantity_factor_and_dimension
from ..dimensions.miscellaneous import dimension_to_si_unit


class Quantity(DimensionSymbol, SymQuantity):  # pylint: disable=too-many-ancestors

    # pylint: disable-next=signature-differs
    def __new__(cls,
        _expr: SupportsFloat = S.One,
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
        return obj

    def __init__(self,
        expr: SupportsFloat = S.One,
        *,
        display_symbol: Optional[str] = None,
        display_latex: Optional[str] = None,
        dimension: Optional[Dimension] = None) -> None:
        (scale, dimension_) = collect_quantity_factor_and_dimension(expr)
        try:
            # if this fails (but it shouldn't), then ``scale`` contains a symbolic sub-expression
            _ = complex(scale)
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

    def split_value_and_unit(self) -> tuple[Expr, Expr]:
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

        return si_value, si_unit

    def _sympystr(self, p: Printer) -> str:
        if "QTY" not in self.display_name:
            return self.display_name

        si_value, si_unit = self.split_value_and_unit()

        return str(p.doprint(si_value * si_unit))

    def _latex(self, p: Printer) -> str:
        if "QTY" not in self.display_latex:
            return self.display_latex

        si_value, si_unit = self.split_value_and_unit()

        return str(p.doprint(si_value * si_unit))


# Allows for some SymPy comparisons, eg Piecewise function
@dispatch(Quantity, Quantity)
def _eval_is_ge(lhs: Quantity, rhs: Quantity) -> bool:
    return scale_factor(lhs) >= scale_factor(rhs)


def subs_list(
    input_: Sequence[SupportsFloat],
    subs_: dict[Expr, SymQuantity],
) -> Sequence[Quantity]:
    return [Quantity(sympify(c, strict=True).subs(subs_)) for c in input_]


def scale_factor(quantity_: SupportsFloat) -> float:
    """
    Extracts the scale factor converted to `float` from the input if it is a quantity. Otherwise
    simply calls `float` on the input.
    """

    if isinstance(quantity_, SymQuantity):
        return float(quantity_.scale_factor)

    return float(quantity_)
