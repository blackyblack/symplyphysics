"""
Relative rocket speed from mass change and effective exhaust speed
==================================================================

The Tsiolkovsky formula determines the speed that an aircraft develops due to the
constant-direction thrust of the rocket engine in the absence of other forces. The generalized
Tsiolkovsky formula is valid for a rocket flying at a speed close to the speed of light.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Tsiolkovsky_rocket_equation#Special_relativity>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)
from symplyphysics.quantities import speed_of_light

speed = symbols.speed
"""
Final :symbols:`speed` of the rocket in the inertial reference frame where the rocket started at
rest.
"""

effective_exhaust_speed = clone_as_symbol(symbols.speed,
    display_symbol="v_e",
    display_latex="v_\\text{e}")
"""
Effective exhaust :symbols:`speed` of the rocket engine.
"""

initial_mass = clone_as_symbol(symbols.mass, subscript="0")
"""
Initial :symbols:`mass` of the rocket.
"""

final_mass = clone_as_symbol(symbols.mass, subscript="1")
"""
Final :symbols:`mass` of the rocket.
"""

law = Eq(final_mass / initial_mass, ((1 - (speed / speed_of_light)) / (1 +
    (speed / speed_of_light)))**(speed_of_light / (2 * effective_exhaust_speed)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(exhaust_velocity_=effective_exhaust_speed,
    initial_mass_=initial_mass,
    final_mass_=final_mass)
@validate_output(speed)
def calculate_speed(exhaust_velocity_: Quantity, initial_mass_: Quantity,
    final_mass_: Quantity) -> Quantity:
    result_expr = solve(law, speed, dict=True)[0][speed]
    result_expr = result_expr.subs({
        effective_exhaust_speed: exhaust_velocity_,
        initial_mass: initial_mass_,
        final_mass: final_mass_
    })
    return Quantity(result_expr)
