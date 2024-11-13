from typing import Self, Any
from sympy import Basic, Symbol as SymSymbol
from ..symbols.symbols import HasDimension, next_name
from ..dimensions import collect_factor_and_dimension


class SingleArgumentOperation(HasDimension, SymSymbol):
    """
    This class is for dimensional symbols that act as an operation with a
    single argument. Its subclasses are intended to be used in custom
    code printers. The argument is placed at position `0` of the `args`
    property.
    """

    def __new__(cls, _argument: Basic, **assumptions: Any) -> Self:
        return SymSymbol.__new__(cls, next_name("SYS"), **assumptions)

    def __init__(self, argument: Basic, **_assumptions: Any) -> None:
        factor, dimension = collect_factor_and_dimension(argument)
        self._args = (factor,)
        HasDimension.__init__(self, dimension)
