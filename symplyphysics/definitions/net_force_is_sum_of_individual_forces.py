"""
Net force is sum of individual forces
=====================================

The net force is the arithmetic sum of forces.

**Conditions:**

#. All force vectors are collinear to each other.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Force#Combining_forces>`__.

#. `Physics LibreTexts, formula 2.2.3 <https://phys.libretexts.org/Courses/University_of_California_Davis/UCD%3A_Physics_9HA__Classical_Mechanics/2%3A_Force/2.2%3A_Effects_of_Force>`__.
"""

from typing import Sequence
from sympy import Eq, Idx, solve
from symplyphysics import (symbols, Quantity, validate_input, validate_output, global_index,
    IndexedSum)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.symbols import clone_as_indexed
from symplyphysics.definitions.vector import superposition_of_forces_is_sum as vector_forces_sum

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, CoordinateVector

net_force = symbols.force
"""
Net :symbols:`force`.
"""

force = clone_as_indexed(symbols.force)
"""
Individual :symbols:`force`.
"""

definition = Eq(net_force, IndexedSum(force[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from the vector form

# Derive the law using 2 forces. Any number of forces can be represented, using 2 of them,
# eg A + B + C = A + (B + C) = Sum(A, Sum(B, C))
_local_index = Idx("local_index_", (1, 2))
_forces_law = definition.subs(global_index, _local_index)
_expected_sum = _forces_law.doit().rhs

# Using one dimensional vectors represents scalar form of the law
_vector_forces = [
    CoordinateVector([force[1], 0, 0], CARTESIAN),
    CoordinateVector([force[2], 0, 0], CARTESIAN),
]
_resultant_vector = vector_forces_sum.law.rhs.subs(global_index, _local_index).doit()
for _index, _force in enumerate(_vector_forces, start=1):
    _resultant_vector = _resultant_vector.subs(vector_forces_sum.force[_index], _force)
_resultant_vector = CoordinateVector.from_expr(_resultant_vector)
assert expr_equals(_resultant_vector.components[0], _expected_sum)
assert _resultant_vector.components[1] == 0
assert _resultant_vector.components[2] == 0


@validate_input(forces_=force)
@validate_output(net_force)
def calculate_resultant_force(forces_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(forces_)))
    forces_law = definition.subs(global_index, local_index)
    forces_law = forces_law.doit()
    solved = solve(forces_law, net_force, dict=True)[0][net_force]
    for i, v in enumerate(forces_):
        solved = solved.subs(force[i + 1], v)
    return Quantity(solved)
