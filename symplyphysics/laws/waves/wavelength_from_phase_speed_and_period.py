"""
Wavelength from phase speed and period
======================================

Wavelength is the spatial period of a periodic wave, i.e. the distance the wave needs to travel
for its shape to repeat. It can also be measured as the distance between two closest points
with the same wave phase.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic import distance_from_constant_velocity as velocity_definition

wavelength = Symbol("wavelength", units.length)
r"""
Wavelength of the wave.

Symbol:
    :code:`lambda`

Latex:
    :math:`\lambda`
"""

phase_velocity = Symbol("phase_velocity", units.velocity)
"""
Phase velocity of the wave.

Symbol:
    :code:`v`
"""

period = Symbol("period", units.time)
"""
Temporal period of the wave.

Symbol:
    :code:`T`
"""

law = Eq(wavelength, phase_velocity * period)
r"""
:code:`lambda = v * T`

Latex:
    .. math::
        \lambda = v T
"""

# Derive the same law from constant velocity motion, assuming wave length is a distance that wave travels
# during oscillation period.

# Prove that derived movement function equals to law.rhs, given initial position = 0
# and phase_velocity is constant_velocity
_constant_velocity_movement_definition = velocity_definition.law.subs({
    velocity_definition.constant_velocity: phase_velocity,
    velocity_definition.movement_time: period,
    velocity_definition.initial_position: 0
})
assert expr_equals(_constant_velocity_movement_definition.rhs, law.rhs)


@validate_input(velocity_=phase_velocity, period_=period)
@validate_output(wavelength)
def calculate_wavelength(velocity_: Quantity, period_: Quantity) -> Quantity:
    applied_definition = solve(law, wavelength, dict=True)[0][wavelength]
    result_expr = applied_definition.subs({
        phase_velocity: velocity_,
        period: period_
    })
    return Quantity(result_expr)
