from typing import Sequence
from sympy import Basic, Expr


# Check if all parameters are being used in trajectory
def contains_parameters(trajectory: Sequence[Expr], parameters: Sequence[Expr]) -> bool:
    free_symbols: set[Basic] = set()
    for component in trajectory:
        free_symbols = free_symbols.union(component.free_symbols)
    for p in parameters:
        if p not in free_symbols:
            return False
    return True


# Parametrized surface should have at least 2 parameters
def is_parametrized_surface(trajectory: Sequence[Expr], parameter1: Expr, parameter2: Expr) -> bool:
    return contains_parameters(trajectory, [parameter1, parameter2])
