from sympy import Expr
from sympy.physics.units import Dimension, Quantity as SymQuantity
from sympy.physics.units.systems.si import SI
from .symbols import DimensionSymbol


class Quantity(DimensionSymbol, SymQuantity):
    def __new__(cls, display_name: str=None, dimension: Dimension=1, scale: Expr=1, **assumptions):
        name = DimensionSymbol.random_name("Q", 8) if display_name is None else DimensionSymbol.random_name(display_name)
        self = super().__new__(cls, name, None, None, None, None, None, False, **assumptions)
        self._dimension = dimension
        self._display_name = name if display_name is None else display_name
        SI.set_quantity_dimension(self, dimension)
        SI.set_quantity_scale_factor(self, scale)
        return self
