from functools import partial
from typing import Any, Optional, Self, Sequence
from sympy import S, Basic, Expr, sympify
from sympy.physics.units import Dimension, Quantity as SymQuantity
from sympy.physics.units.systems.si import SI

from .symbols import DimensionSymbol, next_name
from ..dimensions import collect_factor_and_dimension


class Quantity(DimensionSymbol, SymQuantity):  # pylint: disable=too-many-ancestors

    def __new__(cls,
        _name: Basic,
        _abbrev: Optional[str] = None,
        _latex_repr: Optional[str] = None,
        _pretty_unicode_repr: Optional[str] = None,
        _pretty_ascii_repr: Optional[str] = None,
        _mathml_presentation_repr: Optional[str] = None,
        _is_prefixed: bool = False,
        **assumptions: Any) -> Self:
        name = next_name("QTY")
        obj = SymQuantity.__new__(cls, name, None, None, None, None, None, False, **assumptions)
        return obj

    def __init__(self, expr: Basic | float = S.One, *, dimension: Optional[Dimension] = None):
        (scale, dimension_) = collect_factor_and_dimension(sympify(expr))
        dimension = dimension_ if dimension is None else dimension
        super().__init__(self.name, dimension)
        SI.set_quantity_dimension(self, dimension)
        SI.set_quantity_scale_factor(self, scale)

    # This is required for integration to work properly
    @property
    def func(self) -> partial:
        return partial(Quantity.identity, self)

    def identity(self, *_args: Any) -> Self:
        return self

    def _eval_is_positive(self):
        return self.scale_factor >= 0

    def _eval_Abs(self):
        return self if self.scale_factor >= 0 else (-1 * self)


def list_of_quantities(input_: Sequence[Expr | float], subs_: dict[Expr,
    Quantity]) -> Sequence[Quantity]:
    return [Quantity(sympify(c).subs(subs_)) for c in input_]


def scale_factor(quantity_: Quantity | float) -> float:
    if isinstance(quantity_, Quantity):
        return quantity_.scale_factor
    return quantity_
