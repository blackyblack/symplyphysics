import random
import string
import sympy
from sympy.physics.units import Dimension
from sympy.core.function import UndefinedFunction


class DimensionSymbol:
    _dimension: Dimension = None
    _display_name: str = None

    @staticmethod
    def random_name(base: str=None, digits: int=4) -> str:
        return ("" if base is None else base) + "".join(random.choices(string.digits, k = digits))
    
    @property
    def dimension(self) -> Dimension:
        return self._dimension
    
    @property
    def display_name(self) -> str:
        return self.name if self._display_name is None else self._display_name
    
    def __repr__(self) -> str:
        return self.display_name


class Symbol(DimensionSymbol, sympy.Symbol):    
    def __new__(cls, dimension: Dimension=1, display_name: str=None, **assumptions):
        self = super().__new__(cls, DimensionSymbol.random_name(display_name), **assumptions)
        self._dimension = dimension
        self._display_name = display_name
        return self
    

class Function(DimensionSymbol, UndefinedFunction):    
    def __new__(cls, dimension: Dimension=1, display_name: str=None, **options):
        self = super().__new__(cls, (DimensionSymbol.random_name(display_name)), **options)
        self._dimension = dimension
        self._display_name = display_name
        return self


# Symbol and Function have random name, hence their display is not readable.
# User can provide an expression and a list of it's symbols. All provided
# symbols of this expression will be substituted with their display names.
# The resulting expression should not be reused as it contains very generic
# symbol names that can be modified somewhere else.
def to_printable(expr: sympy.Expr, symbols: list[DimensionSymbol]) -> sympy.Expr:
    new_expr = expr
    for s in symbols:
        symbol = None
        if isinstance(s, sympy.Symbol):
            symbol = sympy.Symbol(s.display_name)
        elif isinstance(s, sympy.FunctionClass):
            symbol = sympy.Function(s.display_name)
        if symbol is None:
            continue
        new_expr = new_expr.subs(s, symbol)
    return new_expr
