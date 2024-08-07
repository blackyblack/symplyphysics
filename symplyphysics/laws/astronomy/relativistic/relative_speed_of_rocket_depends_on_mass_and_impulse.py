"""
Relative speed of rocket depends on mass and impulse
----------------------------------------------------

The Tsiolkovsky formula determines the speed that an aircraft develops under the influence of the thrust
of a rocket engine, unchanged in direction, in the absence of all other forces.
For a rocket flying at a speed close to the speed of light, the generalized Tsiolkovsky formula is valid,
in which the speed of light is present.
"""

from sympy import (Eq, solve)
from sympy.physics.units import speed_of_light
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

speed = Symbol("speed", units.velocity)
"""
Final speed of the rocket

Symbol:
    :code:`V`
"""

exhaust_velocity = Symbol("exhaust_velocity", units.velocity)
"""
Effective exhaust velocity of a rocket engine

Symbol:
    :code:`Ve`

Latex:
    :math:`V_e`
"""

initial_mass = clone_symbol(symbols.basic.mass, "initial_mass")
"""
Initial mass of the rocket

Symbol:
    :code:`M1`

Latex:
    :math:`M_1`
"""

final_mass = clone_symbol(symbols.basic.mass, "final_mass")
"""
Final mass of the rocket

Symbol:
    :code:`M2`

Latex:
    :math:`M_2`
"""

law = Eq(final_mass / initial_mass, ((1 - (speed / speed_of_light)) / (1 +
    (speed / speed_of_light)))**(speed_of_light / (2 * exhaust_velocity)))
r"""
:code:`M2 / M1 = ((1 - V / c) / (1 + V / c))^(c / (2 * Ve))`

Latex:
    .. math::
        \frac{M_2}{M_1} = \left( \frac{1 - \frac{V}{c}}{1 + \frac{V}{c}} \right) ^ {\frac{c}{2 V_e}}
"""


def print_law() -> str:
    return print_expression(law)


@validate_input(exhaust_velocity_=exhaust_velocity,
    initial_mass_=initial_mass,
    final_mass_=final_mass)
@validate_output(speed)
def calculate_speed(exhaust_velocity_: Quantity, initial_mass_: Quantity,
    final_mass_: Quantity) -> Quantity:
    result_expr = solve(law, speed, dict=True)[0][speed]
    result_expr = result_expr.subs({
        exhaust_velocity: exhaust_velocity_,
        initial_mass: initial_mass_,
        final_mass: final_mass_
    })
    return Quantity(result_expr)
