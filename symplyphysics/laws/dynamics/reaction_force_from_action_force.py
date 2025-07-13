"""
Reaction force equals action force
==================================

Newton's third law of motion states that for every action there is an equal and opposite reaction.
Namely, if two bodies exert forces on each other, these forces have the same magnitude but opposite
directions.

**Links:**

#. `Physics LibreTexts, green box <https://phys.libretexts.org/Workbench/PH_245_Textbook_V2/05%3A_Newton's_Laws_of_Motion/5.06%3A_Newtons_Third_Law>`__.

..
    TODO: make vector version of this law
"""

from sympy import (Eq, solve, Idx)
from symplyphysics import (clone_as_symbol, symbols, Quantity, validate_input, validate_output,
    global_index)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.definitions import net_force_is_sum_of_individual_forces as _additive_law
from symplyphysics.laws.conservation import momentum_is_constant as _momentum_conservation_law
from symplyphysics.laws.dynamics import force_is_derivative_of_momentum as _force_momentum_law

action_force = clone_as_symbol(symbols.force, display_symbol="F_12", display_latex="F_{12}")
"""
The projection of the :symbols:`force` that the first body exerts on the second body.
"""

reaction_force = clone_as_symbol(symbols.force, display_symbol="F_21", display_latex="F_{21}")
"""
The projection of the :symbols:`force` that the second body exerts on the first body.
"""

law = Eq(reaction_force, -1 * action_force)
"""
:laws:symbol::

:laws:latex::
"""

# Derive this law from the conservation of momentum and the force-momentum relationship.
# For this law, the closed system consists of the two interacting bodies that only exert forces on
# each other.

_time = _force_momentum_law.time
_momentum = _force_momentum_law.momentum
_force = _force_momentum_law.force

_momentum_conservation_eqn = _momentum_conservation_law.law.subs(
    _momentum_conservation_law.time,
    _time,
).subs(
    _momentum_conservation_law.momentum(_time),
    _momentum(_time),
)

_eqns = (_momentum_conservation_eqn, _force_momentum_law.law)
_net_force_derived = solve(
    _eqns,
    (_momentum(_time).diff(_time), _force(_time)),
    dict=True,
)[0][_force(_time)]

_net_force_expected = _additive_law.definition.rhs.subs(
    global_index,
    Idx("i", (1, 2)),
).doit().subs({
    _additive_law.force[1]: action_force,
    _additive_law.force[2]: reaction_force,
})

_net_force_eqn = Eq(_net_force_expected, _net_force_derived)

_reaction_force_derived = solve(_net_force_eqn, reaction_force)[0]

assert expr_equals(_reaction_force_derived, law.rhs)


@validate_input(force_action_=action_force)
@validate_output(reaction_force)
def calculate_force_reaction(force_action_: Quantity) -> Quantity:
    result_force_expr = solve(law, reaction_force, dict=True)[0][reaction_force]
    result_expr = result_force_expr.subs({action_force: force_action_})
    return Quantity(abs(result_expr))
