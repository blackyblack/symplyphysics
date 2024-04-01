from typing import Sequence
from sympy import Eq, Idx, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    SymbolIndexed,
    SumIndexed,
    global_index,
)

# Description
## For a system of particles, its total rotational inertia is the sum of the rotational
## inertia of each particle of the system.

total_rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
rotational_inertia = SymbolIndexed("rotational_inertia", units.mass * units.length**2)
law = Eq(total_rotational_inertia, SumIndexed(rotational_inertia[global_index], global_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(rotational_inertias_=rotational_inertia)
@validate_output(total_rotational_inertia)
def calculate_rotational_inertia(rotational_inertias_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(rotational_inertias_)))
    rotational_inertia_law = law.subs(global_index, local_index)
    rotational_inertia_law = rotational_inertia_law.doit()
    solved = solve(rotational_inertia_law, total_rotational_inertia,
        dict=True)[0][total_rotational_inertia]
    for i, v in enumerate(rotational_inertias_):
        solved = solved.subs(rotational_inertia[i + 1], v)
    return Quantity(solved)
