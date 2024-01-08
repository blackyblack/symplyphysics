from typing import Sequence
from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, print_expression, Symbol, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## R = sum(F)
## Where:
## F is one of the forces, acting on the body,
## R - resultant (or net) force.

# Conditions:
## - Forces are collinear vectors.

forces = Symbol("forces", units.force)
resultant_force = Symbol("resultant_force", units.force)
definition = Eq(resultant_force, SumArray(forces), evaluate=False)

definition_units_SI = units.force


def print_law() -> str:
    return print_expression(definition)


@validate_input(forces_=forces)
@validate_output(units.force)
def calculate_resultant_force(forces_: Sequence[Quantity]) -> Quantity:
    force_symbols = tuple_of_symbols("force", units.force, len(forces_))
    forces_law = definition.subs(forces, force_symbols).doit()
    solved = solve(forces_law, resultant_force, dict=True)[0][resultant_force]
    for (from_, to_) in zip(force_symbols, forces_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
