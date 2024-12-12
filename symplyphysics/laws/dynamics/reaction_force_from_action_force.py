"""
Reaction force equals action force
==================================

Newton's third law of motion states that for every action there is an equal and opposite reaction.
Namely, if two bodies exert forces on each other, these forces have the same magnitude but opposite
directions.

**Links:**

#. `Physics LibreTexts, green box <https://phys.libretexts.org/Workbench/PH_245_Textbook_V2/05%3A_Newton's_Laws_of_Motion/5.06%3A_Newtons_Third_Law>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (clone_as_symbol, symbols, Quantity, validate_input, validate_output)

action_force = clone_as_symbol(symbols.force, display_symbol="F_12", display_latex="F_{12}")
"""
The projection of the :symbols:`force` that the first body exerts on the second body.
"""

reaction_force = clone_as_symbol(symbols.force, display_symbol="F_21", display_latex="F_{21}")
"""
The projection of the :symbols:`force` that the second body exerts on the first body.
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
