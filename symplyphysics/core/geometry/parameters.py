from typing import Sequence, TypeAlias
from sympy import Basic, Expr
from symplyphysics.core.dimensions import ScalarValue

ParameterType: TypeAlias = Expr
ParameterLimits: TypeAlias = tuple[ParameterType, ScalarValue, ScalarValue]


# Check if all parameters are being used in trajectory
def contains_parameters(trajectory: Sequence[Expr], parameters: Sequence[ParameterType]) -> bool:
    free_symbols: set[Basic] = set()
    for component in trajectory:
        free_symbols = free_symbols.union(component.free_symbols)
    for p in parameters:
        if p not in free_symbols:
            return False
    return True


# Parametrized surface should have at least 2 parameters
def is_parametrized_surface(trajectory: Sequence[Expr], parameter1: ParameterType,
    parameter2: ParameterType) -> bool:
    return contains_parameters(trajectory, [parameter1, parameter2])
