from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## For a system of particles, its total rotational inertia is the sum of the rotational
## inertia of each particle of the system.

rotational_inertias = Symbol("rotational_inertias", units.mass * units.length**2)
total_rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)

law = Eq(total_rotational_inertia, SumArray(rotational_inertias), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(rotational_inertias_=rotational_inertias)
@validate_output(total_rotational_inertia)
def calculate_rotational_inertia(rotational_inertias_: list[Quantity]) -> Quantity:
    rotational_inertia_symbols = tuple_of_symbols("rotational_inertia",
        units.mass * units.length**2, len(rotational_inertias_))
    rotational_inertia_law = law.subs(rotational_inertias, rotational_inertia_symbols).doit()
    solved = solve(rotational_inertia_law, total_rotational_inertia)[0]
    for (symbol, value) in zip(rotational_inertia_symbols, rotational_inertias_):
        solved = solved.subs(symbol, value)
    return Quantity(solved)
