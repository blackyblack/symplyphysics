"""
Friction force from normal force
================================

The *friction* force is tangential interaction between two objects, which impedes their relative movement.
It is proportional to the normal force between the two objects.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

friction_force = clone_as_symbol(symbols.force,
    display_symbol="F_fr",
    display_latex="F_\\text{fr}")
"""
The friction :symbols:`force`.
"""

coefficient_of_friction = symbols.coefficient_of_friction
r"""
The :symbols:`coefficient_of_friction` between the two objects.
"""

normal_force = clone_as_symbol(
    symbols.force,
    display_symbol="N",
    display_latex="N",
)
"""
The normal reaction :symbols:`force` from one object to another.
"""

law = Eq(friction_force, coefficient_of_friction * normal_force)
r"""
:code:`F_fr = mu * N`

Latex:
    .. math::
        F_\text{fr} = \mu N
"""


@validate_input(coefficient_of_friction_=coefficient_of_friction, normal_reaction_=normal_force)
@validate_output(friction_force)
def calculate_friction_force(coefficient_of_friction_: float,
    normal_reaction_: Quantity) -> Quantity:
    result_expr = solve(law, friction_force, dict=True)[0][friction_force]
    friction_force_applied = result_expr.subs({
        coefficient_of_friction: coefficient_of_friction_,
        normal_force: normal_reaction_
    })
    return Quantity(friction_force_applied)
