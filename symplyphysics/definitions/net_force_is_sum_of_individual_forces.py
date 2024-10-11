"""
Net force is sum of individual forces
=====================================

The net force is the arithmetic sum of forces.

**Conditions:**

#. All force vectors are collinear to each other.
"""

from typing import Sequence
from sympy import (Eq, Idx, solve)
from symplyphysics import (symbols, units, Quantity, validate_input, validate_output, global_index,
    SumIndexed, Vector)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.symbols import clone_as_indexed
from symplyphysics.definitions.vector import superposition_of_forces_is_sum as vector_forces_sum

net_force = symbols.force
"""
Net :symbols:`force`.
"""

force = clone_as_indexed(symbols.force, display_symbol="F[i]", display_latex="F_i")
"""
Individual :symbols:`force`.
"""

definition = Eq(net_force, SumIndexed(force[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from the vector form

# Derive the law using 2 forces. Any number of forces can be represented, using 2 of them,
# eg A + B + C = A + (B + C) = Sum(A, Sum(B, C))
_local_index_ = Idx("local_index_", (1, 2))
_forces_law_ = definition.subs(global_index, _local_index_)
_expected_sum = _forces_law_.doit().rhs

# Using one dimensional vectors represents scalar form of the law
_vector_forces = [Vector([force[1]]), Vector([force[2]])]
_resultant_vector = vector_forces_sum.superposition_law(_vector_forces)
assert expr_equals(_resultant_vector.components[0], _expected_sum)
for component in _resultant_vector.components[1:]:
    assert component == 0


@validate_input(forces_=force)
@validate_output(units.force)
def calculate_resultant_force(forces_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(forces_)))
    forces_law = definition.subs(global_index, local_index)
    forces_law = forces_law.doit()
    solved = solve(forces_law, net_force, dict=True)[0][net_force]
    for i, v in enumerate(forces_):
        solved = solved.subs(force[i + 1], v)
    return Quantity(solved)
