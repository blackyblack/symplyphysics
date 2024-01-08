from typing import Sequence
from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, print_expression, Symbol, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols
from symplyphysics.core.vectors.vectors import Vector
from symplyphysics.definitions.vector import superposition_of_forces_is_sum as vector_forces_sum

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

# Derive the same law from the vector form

# Derive the law using 2 forces. Any number of forces can be represented, using 2 of them,
# eg A + B + C = A + (B + C) = Sum(A, Sum(B, C))
force_symbols_ = tuple_of_symbols("force", units.force, 2)
(force1, force2) = force_symbols_
expected_sum = definition.subs(forces, force_symbols_).doit().rhs

# Using one dimensional vectors represents scalar form of the law
vector_forces = [Vector([force1]), Vector([force2])]
resultant_vector = vector_forces_sum.superposition_law(vector_forces)
assert len(resultant_vector.components) == 1
assert expr_equals(resultant_vector.components[0], expected_sum)


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
