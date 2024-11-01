"""
Wavelength from phase speed and period
======================================

Wavelength is the spatial period of a periodic wave, i.e. the distance the wave needs to travel
for its shape to repeat. It can also be measured as the distance between two closest points
with the same wave phase.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics import (
    position_via_constant_speed_and_time as velocity_definition,
)

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the wave.
"""

phase_velocity = symbols.phase_speed
"""
:symbols:`phase_speed` of the wave.
"""

period = symbols.period
"""
Temporal :symbols:`period` of the wave.
"""

law = Eq(wavelength, phase_velocity * period)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from constant velocity motion, assuming wave length is a distance that wave travels
# during oscillation period.

# Prove that derived movement function equals to law.rhs, given initial position = 0
# and phase_velocity is constant_velocity
_constant_velocity_movement_definition = velocity_definition.law.subs({
    velocity_definition.speed: phase_velocity,
    velocity_definition.time: period,
    velocity_definition.initial_position: 0
})
assert expr_equals(_constant_velocity_movement_definition.rhs, law.rhs)


@validate_input(velocity_=phase_velocity, period_=period)
@validate_output(wavelength)
def calculate_wavelength(velocity_: Quantity, period_: Quantity) -> Quantity:
    applied_definition = solve(law, wavelength, dict=True)[0][wavelength]
    result_expr = applied_definition.subs({phase_velocity: velocity_, period: period_})
    return Quantity(result_expr)
