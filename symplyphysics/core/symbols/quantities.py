from functools import partial
from sympy import S, Expr
from sympy.physics.units import Dimension, Quantity as SymQuantity
from sympy.physics.units.systems.si import SI
from ..expr_to_quantity import collect_factor_and_dimension
from .symbols import DimensionSymbol


class Quantity(DimensionSymbol, SymQuantity):
    def __new__(cls, scale: Expr=1, *, display_name: str=None, dimension: Dimension=None, **assumptions):
        (_, dimension_) = collect_factor_and_dimension(scale)
        dimension = dimension_ if dimension is None else dimension
        name = DimensionSymbol.random_name("Q", 8) if display_name is None else DimensionSymbol.random_name(display_name)
        self = super().__new__(cls, name, None, None, None, None, None, False, **assumptions)
        self._dimension = dimension
        self._display_name = name if display_name is None else display_name
        SI.set_quantity_dimension(self, dimension)
        SI.set_quantity_scale_factor(self, scale)
        return self
    
    # This is required for integration to work properly
    @property
    def func(self):
        return partial(Quantity.identity, self)

    def identity(self, *_args):
        return self
    

Dimensionless = Dimension(S.One)
