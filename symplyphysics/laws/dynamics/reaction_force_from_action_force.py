"""
Reaction force equals action force
==================================

Newton's third law of motion states that for every action there is an equal and opposite reaction.
Namely, if two bodies exert forces on each other, these forces have the same magnitude but opposite
directions.
"""

from sympy import (Eq, solve)
from symplyphysics import (clone_symbol, symbols, Quantity, validate_input, validate_output)

action_force = clone_symbol(symbols.dynamics.force, "action_force")
r"""
The projection of the :attr:`~symplyphysics.symbols.dynamics.force` that the first body exerts on the second body.

Symbol:
    :code:`F_12`

Latex:
    :math:`F_{12}`
"""

reaction_force = clone_symbol(symbols.dynamics.force, "reaction_force")
r"""
The projection of the :attr:`~symplyphysics.symbols.dynamics.force` that the second body exerts on the first body.

Symbol:
    :code:`F_21`

Latex:
    :math:`F_{21}`
"""

law = Eq(reaction_force, -1 * action_force)
r"""
:code:`F_12 = -1 * F_21`

Latex:
    .. math::
        F_{12} = - F_{21}
"""


@validate_input(force_action_=action_force)
@validate_output(reaction_force)
def calculate_force_reaction(force_action_: Quantity) -> Quantity:
    result_force_expr = solve(law, reaction_force, dict=True)[0][reaction_force]
    result_expr = result_force_expr.subs({action_force: force_action_})
    return Quantity(abs(result_expr))
