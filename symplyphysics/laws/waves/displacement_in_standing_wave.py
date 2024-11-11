r"""
Displacement in standing wave
=============================

A *standing*, or *stationary*, *wave* is the result of the interference of two identical waves
moving in opposite directions.

**Notes:**

#. In this law we assume the standing wave to be composed of two identical traveling sinusoidal
   waves of the form :math:`u_\text{max} \sin(k x \pm \omega t)`
#. A standing wave is no longer a traveling one because it doesn't move in a single direction.
"""

from sympy import Eq, sin, cos, symbols as sym_symbols
from symplyphysics import Quantity, validate_input, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.quantity_decorator import validate_output_same
from symplyphysics.laws.waves import phase_of_traveling_wave as phase_law

total_displacement = sym_symbols("u")
"""
Displacement of the resulting wave.

Symbol:
    :code:`u`

Latex:
    :math:`u`
"""

amplitude = sym_symbols("u_max")
r"""
Amplitude of the interfering waves.

Symbol:
    :code:`u_max`

Latex:
    :math:`u_\text{max}`
"""

angular_wavenumber = symbols.angular_wavenumber
"""
:symbols:`angular_wavenumber` of the interfering waves.
"""

position = symbols.position
"""
:symbols:`position`, or spatial coordinate.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the interfering waves.
"""

time = symbols.time
"""
:symbols:`time`.
"""

law = Eq(total_displacement,
    (2 * amplitude) * sin(angular_wavenumber * position) * cos(angular_frequency * time))
r"""
:code:`u = 2 * u_max * sin(k * x) * cos(w * t)`

Latex:
    .. math::
        u = 2 u_\text{max} \sin(k x) \cos(\omega t)
"""

# Derive from the sum of two traveling waves

_wave_phase = phase_law.law.rhs.subs({
    phase_law.angular_wavenumber: angular_wavenumber,
    phase_law.position: position,
    phase_law.angular_frequency: angular_frequency,
    phase_law.time: time,
})

# The form of the waves is taken from the note above
_forward_wave = amplitude * sin(_wave_phase)

_backward_wave = _forward_wave.subs(angular_frequency, -1 * angular_frequency)

_sum_of_waves = (_forward_wave + _backward_wave).simplify()

assert expr_equals(_sum_of_waves, law.rhs)


@validate_input(
    angular_wavenumber_=angular_wavenumber,
    position_=position,
    angular_frequency_=angular_frequency,
    time_=time,
)
@validate_output_same("amplitude_")
def calculate_displacement(
    amplitude_: Quantity,
    angular_wavenumber_: Quantity,
    position_: Quantity,
    angular_frequency_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        amplitude: amplitude_,
        angular_wavenumber: angular_wavenumber_,
        position: position_,
        angular_frequency: angular_frequency_,
        time: time_,
    })
    return Quantity(result)
