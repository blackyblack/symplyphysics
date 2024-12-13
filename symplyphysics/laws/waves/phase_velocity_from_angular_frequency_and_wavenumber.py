"""
Phase speed from angular frequency and wavenumber
=================================================

*Phase speed* of a wave is the rate at which the wave propagates in a medium.
It is the speed at which the phase of one frequency component of the wave travels.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Phase_velocity>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import angular_wavenumber_is_inverse_wavelength as wavenumber_def
from symplyphysics.laws.waves import wavelength_from_phase_speed_and_period as velocity_law
from symplyphysics.definitions import period_from_angular_frequency as frequency_law

phase_speed = symbols.phase_speed
"""
:symbols:`phase_speed` of the wave.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the wave.
"""

angular_wavenumber = symbols.angular_wavenumber
"""
:symbols:`angular_wavenumber` of the wave.
"""

law = Eq(phase_speed, angular_frequency / angular_wavenumber)
"""
:laws:symbol::

:laws:latex::
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
