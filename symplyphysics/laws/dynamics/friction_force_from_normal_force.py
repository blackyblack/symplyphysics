"""
Friction force from normal force
================================

The *friction* force is tangential interaction between two objects, which impedes their relative movement.
It is proportional to the normal force between the two objects.

**Links**

#. `Physics LibreTexts <https://phys.libretexts.org/Courses/Tuskegee_University/Algebra_Based_Physics_I/04%3A_Dynamics-_Force_and_Newton%27s_Laws_of_Motion/4.07%3A_Friction>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

friction_force = clone_as_symbol(symbols.force, display_symbol="F_fr", display_latex="F_\\text{fr}")
"""
The friction :symbols:`force`.
"""

friction_coefficient = symbols.friction_coefficient
r"""
The :symbols:`friction_coefficient` between the two objects.
"""

normal_force = clone_as_symbol(
    symbols.force,
    display_symbol="N",
    display_latex="N",
)
"""
The normal reaction :symbols:`force` from one object to another.
"""

law = Eq(friction_force, friction_coefficient * normal_force)
r"""
:code:`F_fr = mu * N`

Latex:
    .. math::
        F_\text{fr} = \mu N
"""


@validate_input(friction_coefficient_=friction_coefficient, normal_reaction_=normal_force)
@validate_output(friction_force)
def calculate_friction_force(friction_coefficient_: float, normal_reaction_: Quantity) -> Quantity:
    result_expr = solve(law, friction_force, dict=True)[0][friction_force]
    friction_force_applied = result_expr.subs({
        friction_coefficient: friction_coefficient_,
        normal_force: normal_reaction_
    })
    return Quantity(friction_force_applied)
