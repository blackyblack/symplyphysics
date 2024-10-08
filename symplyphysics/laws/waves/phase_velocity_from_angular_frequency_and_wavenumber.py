"""
Phase speed from angular frequency and wavenumber
=================================================

*Phase speed* of a wave is the rate at which the wave propagates in a medium.
It is the speed at which the phase of one frequency component of the wave travels.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    angle_type,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import angular_wavenumber_is_inverse_wavelength as wavenumber_def
from symplyphysics.laws.waves import wavelength_from_phase_speed_and_period as velocity_law
from symplyphysics.definitions import period_from_angular_frequency as frequency_law

phase_speed = Symbol("phase_speed", units.length / units.time, real=True)
"""
Phase speed of the wave.

Symbol:
    :code:`v`
"""

angular_frequency = Symbol("angular_frequency", angle_type / units.time, positive=True)
r"""
Angular frequency of the wave.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

angular_wavenumber = Symbol("angular_wavenumber", angle_type / units.length, positive=True)
"""
Angular wavenumber of the wave.

Symbol:
    :code:`k`
"""

law = Eq(phase_speed, angular_frequency / angular_wavenumber)
r"""
:code:`v = w / k`

Latex:
    .. math::
        v = \frac{\omega}{k}
"""

# Derive from definition of phase velocity via wavelength and time period

_wavelength = solve(wavenumber_def.definition,
    wavenumber_def.wavelength)[0].subs(wavenumber_def.angular_wavenumber, angular_wavenumber)

_time_period = solve(frequency_law.law,
    frequency_law.period)[0].subs(frequency_law.angular_frequency, angular_frequency)

_phase_velocity = solve(velocity_law.law, velocity_law.phase_velocity)[0].subs({
    velocity_law.wavelength: _wavelength,
    velocity_law.period: _time_period,
})

assert expr_equals(law.rhs, _phase_velocity)


@validate_input(
    angular_frequency_=angular_frequency,
    angular_wavenumber_=angular_wavenumber,
)
@validate_output(phase_speed)
def calculate_phase_velocity(
    angular_frequency_: Quantity,
    angular_wavenumber_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        angular_frequency: angular_frequency_,
        angular_wavenumber: angular_wavenumber_,
    })
    return Quantity(result)
