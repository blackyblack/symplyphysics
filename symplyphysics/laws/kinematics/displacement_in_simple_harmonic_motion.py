"""
Displacement in simple harmonic motion
======================================

Any motion that repeats at regular intervals is called periodic, or harmonic, motion.
Simple harmonic motion is a particular type of repeated motion that is a sinusoidal
function of time.

**Note:**

#. This law is also applicable for any physical quantity that changes its value in
   a repeating harmonic manner, therefore :symbols:`any_quantity` is used.
"""

from sympy import Eq, cos, dsolve
from symplyphysics import (
    Quantity,
    validate_input,
    symbols,
    clone_as_symbol,
    clone_as_function,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.quantity_decorator import validate_output_same
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.definitions import harmonic_oscillator_is_second_derivative_equation as harmonic_eqn

displacement = clone_as_symbol(
    symbols.any_quantity,
    display_symbol="q",
    display_latex="q",
)
"""
Displacement from rest, usually a function of :symbols:`time`. See :symbols:`any_quantity`.
"""

amplitude = clone_as_symbol(
    symbols.any_quantity,
    display_symbol="q_max",
    display_latex="q_\\text{max}",
)
"""
The maximum absolute value of the :attr:`~displacement`. See :symbols:`any_quantity`.
"""

angular_frequency = clone_as_symbol(symbols.angular_frequency, positive=True)
"""
:symbols:`angular_frequency` of oscillations.
"""

phase_shift = symbols.phase_shift
"""
:symbols:`phase_shift` of oscillations, which is the phase at :math:`t = 0`.
"""

time = symbols.time
"""
:symbols:`time` at which :attr:`~displacement` is measured.
"""

law = Eq(displacement, amplitude * cos(angular_frequency * time + phase_shift))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from [oscillator equation](../../definitions/harmonic_oscillator_is_second_derivative_equation.py)

_displacement = clone_as_function(displacement, [time])

_eqn = harmonic_eqn.definition.replace(harmonic_eqn.displacement, _displacement).subs({
    harmonic_eqn.time: time,
    harmonic_eqn.angular_frequency: angular_frequency,
})

_initial_position = law.rhs.subs(time, 0)
_initial_velocity = law.rhs.diff(time).subs(time, 0)

_dsolved = dsolve(_eqn,
    _displacement(time),
    ics={
    _displacement(0): _initial_position,
    _displacement(time).diff(time).subs(time, 0): _initial_velocity,
    }).rhs

assert expr_equals(law.rhs, _dsolved)


@validate_input(
    angular_frequency_=angular_frequency,
    phase_lag_=phase_shift,
    time_=time,
)
@validate_output_same("amplitude_")
def calculate_displacement(
    amplitude_: Quantity,
    angular_frequency_: Quantity,
    phase_lag_: Quantity | float,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        amplitude: amplitude_,
        angular_frequency: angular_frequency_,
        phase_shift: scale_factor(phase_lag_),
        time: time_,
    })
    return Quantity(result)
