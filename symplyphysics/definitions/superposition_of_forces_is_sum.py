from typing import Sequence
from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, print_expression, Symbol, validate_input,
    validate_output, SymbolIndexed, global_index, SumIndexed, Vector)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions.vector import superposition_of_forces_is_sum as vector_forces_sum

# Description
## R = sum(F)
## Where:
## F is one of the forces, acting on the body,
## R - resultant (or net) force.

# Conditions:
## - Forces are collinear vectors.

resultant_force = Symbol("resultant_force", units.force)
force = SymbolIndexed("force", units.force)
definition = Eq(resultant_force, SumIndexed(force[global_index], global_index))

definition_units_SI = units.force

# Derive the same law from the vector form

# Derive the law using 2 forces. Any number of forces can be represented, using 2 of them,
# eg A + B + C = A + (B + C) = Sum(A, Sum(B, C))
local_index_ = Idx("local_index_", (1, 2))
forces_law_ = definition.subs(global_index, local_index_)
expected_sum = forces_law_.doit().rhs

# Using one dimensional vectors represents scalar form of the law
vector_forces = [Vector([force[1]]), Vector([force[2]])]
resultant_vector = vector_forces_sum.superposition_law(vector_forces)
assert len(resultant_vector.components) == 3
assert expr_equals(resultant_vector.components[0], expected_sum)
assert resultant_vector.components[1] == 0
assert resultant_vector.components[2] == 0


def print_law() -> str:
    return print_expression(definition)


@validate_input(forces_=force)
@validate_output(units.force)
def calculate_resultant_force(forces_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(forces_)))
    forces_law = definition.subs(global_index, local_index)
    forces_law = forces_law.doit()
    solved = solve(forces_law, resultant_force, dict=True)[0][resultant_force]
    for i, v in enumerate(forces_):
        solved = solved.subs(force[i + 1], v)
    return Quantity(solved)
