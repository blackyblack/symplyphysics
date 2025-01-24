"""
Stiffness of two parallel springs
=================================

If two springs are side-to-side to one another, i.e. connected in parallel, the total
stiffness of the system of springs is equal to the sum of the stiffnesses of each
spring.

**Condition:**

#. Springs must be Hookean, or linear-response, i.e. obey the Hooke's law.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Series_and_parallel_springs#Formulas>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    Vector,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions.vector import superposition_of_forces_is_sum as superposition_law
from symplyphysics.laws.dynamics.springs import spring_reaction_is_proportional_to_deformation as hookes_law

total_stiffness = symbols.stiffness
"""
Total :symbols:`stiffness`.
"""

first_stiffness = clone_as_symbol(symbols.stiffness, subscript="1")
"""
:symbols:`stiffness` of the first spring.
"""

second_stiffness = clone_as_symbol(symbols.stiffness, subscript="2")
"""
:symbols:`stiffness` of the second spring.
"""

law = Eq(total_stiffness, first_stiffness + second_stiffness)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from Hooke's law
## In the case of parallel connection, the forces acting on both springs add up while
## the springs deform equally.

_deformation = symbols.deformation

_force_expr = solve(hookes_law.law,
    hookes_law.spring_reaction)[0].subs(hookes_law.deformation, _deformation)

_first_force = _force_expr.subs(hookes_law.stiffness, first_stiffness)
_second_force = _force_expr.subs(hookes_law.stiffness, second_stiffness)

_total_force_hooke = _force_expr.subs(hookes_law.stiffness, total_stiffness)

_first_force_vector = Vector([_first_force])
_second_force_vector = Vector([_second_force])
_total_force_vector = superposition_law.superposition_law(
    [_first_force_vector, _second_force_vector])
for component in _total_force_vector.components[1:]:
    assert expr_equals(component, 0)
_total_force_added = _total_force_vector.components[0]

_total_stiffness_derived = solve(
    Eq(_total_force_hooke, _total_force_added),
    total_stiffness,
)[0]

assert expr_equals(_total_stiffness_derived, law.rhs)


@validate_input(
    first_stiffness_=first_stiffness,
    second_stiffness_=second_stiffness,
)
@validate_output(total_stiffness)
def calculate_total_stiffness(
    first_stiffness_: Quantity,
    second_stiffness_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        first_stiffness: first_stiffness_,
        second_stiffness: second_stiffness_,
    })
    return Quantity(result)
