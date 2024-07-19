"""
Friction force from normal force
================================

The *friction* force is tangential interaction between two objects, which impedes there relative movement.
It is proportional to the normal force between the two objects.
"""

from sympy import (Eq, solve)
from symplyphysics import (clone_symbol, symbols, Quantity, Symbol, dimensionless,
    validate_input, validate_output)

friction_force = clone_symbol(symbols.dynamics.force, "friction_force")
r"""
The friction :attr:`~symplyphysics.symbols.dynamics.force`.

Symbol:
    F_friction

Latex:
    :math:`F_\text{fr}`
"""

friction_factor = Symbol("friction_factor", dimensionless)
r"""
The friction factor between the two objects.

Symbol:
    mu

Latex:
    :math:`\mu`
"""

normal_reaction = clone_symbol(symbols.dynamics.force, "normal_reaction")
"""
The normal reaction :attr:`~symplyphysics.symbols.dynamics.force` from one object to another.

Symbol:
    N
"""

law = Eq(friction_force, friction_factor * normal_reaction)
r"""
F_friction = mu * N

Latex:
    :math:`F_\text{fr} = \mu N`
"""


@validate_input(friction_factor_=friction_factor, normal_reaction_=normal_reaction)
@validate_output(friction_force)
def calculate_friction_force(friction_factor_: float, normal_reaction_: Quantity) -> Quantity:
    result_expr = solve(law, friction_force, dict=True)[0][friction_force]
    friction_force_applied = result_expr.subs({
        friction_factor: friction_factor_,
        normal_reaction: normal_reaction_
    })
    return Quantity(friction_force_applied)
