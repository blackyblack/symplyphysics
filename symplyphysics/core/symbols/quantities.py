from functools import partial
from typing import Any, Optional, Self, Sequence
from sympy import S, Basic, Expr, sympify, Abs
from sympy.physics.units import Dimension, Quantity as SymQuantity
from sympy.physics.units.systems.si import SI
from sympy.multipledispatch import dispatch

from .symbols import DimensionSymbolNew, next_name
from ..dimensions import collect_factor_and_dimension
from ..errors import UnitsError


class Quantity(DimensionSymbolNew, SymQuantity):  # pylint: disable=too-many-ancestors

    # pylint: disable-next=signature-differs
    def __new__(cls,
        _expr: Basic | float = S.One,
        *,
        display_symbol: Optional[str] = None,
        display_latex: Optional[str] = None,
        **assumptions: Any) -> Self:
        name = next_name("QTY")
        # Latex symbol is set in SymPy Quantity, not in DimensionSymbol, due to Latex printer
        # specifics
        display_symbol = name if display_symbol is None else display_symbol
        display_latex = display_symbol if display_latex is None else display_latex
        obj = SymQuantity.__new__(cls, name, None, display_latex, None, None, None, False,
            **assumptions)
        return obj

    def __init__(self,
        expr: Basic | float = S.One,
        *,
        display_symbol: Optional[str] = None,
        display_latex: Optional[str] = None,
        dimension: Optional[Dimension] = None) -> None:
        (scale, dimension_) = collect_factor_and_dimension(sympify(expr))
        if scale.free_symbols:
            raise UnitsError(f"Argument '{expr}' to function 'Quantity()' should "
                f"not contain free symbols")
        dimension = dimension_ if dimension is None else dimension
        display_symbol = str(self.name) if display_symbol is None else display_symbol
        super().__init__(display_symbol, dimension, display_latex=display_latex)
        SI.set_quantity_dimension(self, dimension)
        SI.set_quantity_scale_factor(self, scale)

    # This is required for integration to work properly
    @property
    def func(self) -> partial:
        return partial(Quantity.identity, self)

    def identity(self, *_args: Any) -> Self:
        return self

    def _eval_is_positive(self) -> bool:
        return self.scale_factor >= 0

    def _eval_Abs(self) -> Self:
        return self.__class__(Abs(self.scale_factor), dimension=self.dimension)


# Allows for some SymPy comparisons, eg Piecewise function
@dispatch(Quantity, Quantity)
def _eval_is_ge(lhs: Quantity, rhs: Quantity) -> bool:
    return lhs.scale_factor >= rhs.scale_factor


def subs_list(input_: Sequence[Expr | float], subs_: dict[Expr, Quantity]) -> Sequence[Quantity]:
    return [Quantity(sympify(c).subs(subs_)) for c in input_]


def scale_factor(quantity_: Quantity | float) -> float:
    return quantity_.scale_factor if isinstance(quantity_, Quantity) else quantity_


def evaluate_quantity(quantity: Expr, **kwargs: Any) -> Quantity:
    if not isinstance(quantity, Quantity):
        quantity = Quantity(quantity)
    scale_factor_ = quantity.scale_factor.evalf(**kwargs)
    dimension = quantity.dimension
    return Quantity(scale_factor_, dimension=dimension)
