"""
Friction force from normal force
================================

The *friction* force is tangential interaction between two objects, which impedes their relative movement.
It is proportional to the normal force between the two objects.
"""

from sympy import (Eq, solve)
from symplyphysics import (clone_symbol, symbols, Quantity, Symbol, dimensionless,
    validate_input, validate_output)

friction_force = clone_symbol(symbols.dynamics.force, "friction_force")
r"""
The friction :attr:`~symplyphysics.symbols.dynamics.force`.

Symbol:
    :code:`F_fr`

Latex:
    :math:`F_\text{fr}`
"""

coefficient_of_friction = Symbol("coefficient_of_friction", dimensionless)
r"""
The coefficient of friction between the two objects.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

normal_force = clone_symbol(symbols.dynamics.force, "normal_force")
"""
The normal reaction :attr:`~symplyphysics.symbols.dynamics.force` from one object to another.

Symbol:
    :code:`N`
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
def calculate_friction_force(coefficient_of_friction_: float, normal_reaction_: Quantity) -> Quantity:
    result_expr = solve(law, friction_force, dict=True)[0][friction_force]
    friction_force_applied = result_expr.subs({
        coefficient_of_friction: coefficient_of_friction_,
        normal_force: normal_reaction_
    })
    return Quantity(friction_force_applied)
