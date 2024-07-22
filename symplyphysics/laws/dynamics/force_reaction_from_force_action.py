"""
Force reaction from force action
================================

Newton's third law of motion states that for every action there is an equal and opposite reaction.
Namely, if two bodies exert forces on each other, these forces have the same magnitude but opposite
directions.
"""

from sympy import (Eq, solve)
from symplyphysics import (clone_symbol, symbols, Quantity, validate_input,
    validate_output)

force_action = clone_symbol(symbols.dynamics.force, "force_action")
r"""
The projection of the :attr:`~symplyphysics.symbols.dynamics.force` that the first body exerts on the second body.

Symbol:
    F_12

Latex:
    :math:`F_{12}`
"""

force_reaction = clone_symbol(symbols.dynamics.force, "force_reaction")
r"""
The projection of the :attr:`~symplyphysics.symbols.dynamics.force` that the second body exerts on the first body.

Symbol:
    F_21

Latex:
    :math:`F_{21}`
"""

law = Eq(force_reaction, -1 * force_action)
r"""
F_12 = -1 * F_21

Latex:
    .. math::
        F_{12} = - F_{21}
"""


@validate_input(force_action_=force_action)
@validate_output(force_reaction)
def calculate_force_reaction(force_action_: Quantity) -> Quantity:
    result_force_expr = solve(law, force_reaction, dict=True)[0][force_reaction]
    result_expr = result_force_expr.subs({force_action: force_action_})
    return Quantity(abs(result_expr))
