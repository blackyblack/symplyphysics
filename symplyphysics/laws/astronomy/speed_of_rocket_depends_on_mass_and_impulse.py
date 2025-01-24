"""
Rocket speed from mass and impulse
==================================

The Tsiolkovsky formula determines the speed that an aircraft develops under the influence of the thrust
of a rocket engine, unchanged in direction, in the absence of all other forces.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Tsiolkovsky_rocket_equation>`__.
"""

from sympy import (Eq, solve, log)
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

speed_change = clone_as_symbol(symbols.speed, display_symbol="Delta(v)", display_latex="\\Delta v")
"""
Maximum change in :symbols:`speed` of the rocket.
"""

effective_exhaust_speed = clone_as_symbol(symbols.speed, display_symbol="v_e", display_latex="v_\\text{e}")
"""
Effective exhaust :symbols:`speed`, or specific impulse (i.e. impulse per unit mass). See
`this Wikipedia paragraph <https://en.wikipedia.org/wiki/Specific_impulse#Specific_impulse_as_effective_exhaust_velocity>`__
for more information.
"""

initial_mass = clone_as_symbol(symbols.mass, subscript="0")
"""
Initial :symbols:`mass` of the rocket.
"""

final_mass = clone_as_symbol(symbols.mass, subscript="1")
"""
Final :symbols:`mass` of the rocket.
"""

law = Eq(speed_change, effective_exhaust_speed * log(initial_mass / final_mass))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(specific_impulse_=effective_exhaust_speed,
    initial_mass_=initial_mass,
    final_mass_=final_mass)
@validate_output(speed_change)
def calculate_speed(specific_impulse_: Quantity, initial_mass_: Quantity,
    final_mass_: Quantity) -> Quantity:
    result_expr = solve(law, speed_change, dict=True)[0][speed_change]
    result_expr = result_expr.subs({
        effective_exhaust_speed: specific_impulse_,
        initial_mass: initial_mass_,
        final_mass: final_mass_
    })
    return Quantity(result_expr)
