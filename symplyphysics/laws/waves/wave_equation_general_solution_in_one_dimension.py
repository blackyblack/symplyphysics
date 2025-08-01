"""
General solution to wave equation in one dimension
==================================================

Any function of a single variable can be a solution of the wave equation if it depends
on the phase of the wave, i.e. the position and time variables only appear in its
expression in the form of the :doc:`wave phase <laws.waves.phase_of_traveling_wave>`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Wave_equation#General_solution>`__.
"""

from sympy import Eq, cos
from symplyphysics import Quantity, validate_input, symbols, Symbol, Function
from symplyphysics.core.dimensions import any_dimension
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.quantity_decorator import validate_output_same
from symplyphysics.core.symbols.quantities import scale_factor

from symplyphysics.definitions import wave_equation_in_one_dimension as wave_eqn
from symplyphysics.laws.waves import (
    phase_of_traveling_wave as phase_law,
    phase_velocity_from_angular_frequency_and_wavenumber as phase_velocity_law,
)

displacement = Symbol("u", any_dimension)
"""
Displacement from rest in the wave. Usually depends on position and time.
"""

wave_phase = symbols.phase
"""
:symbols:`phase` of the wave.
"""

solution = Function("f", (wave_phase,), dimension=any_dimension)
"""
One-argument solution function of the wave equation.
"""

law = Eq(displacement, solution(wave_phase))
"""
:laws:symbol::

:laws:latex::
"""

# Prove this is a solution of the wave equation

_angular_frequency = symbols.angular_frequency
_angular_wavenumber = symbols.angular_wavenumber
_position = symbols.position
_time = symbols.time

_phase_velocity = phase_velocity_law.law.rhs.subs({
    phase_velocity_law.angular_frequency: _angular_frequency,
    phase_velocity_law.angular_wavenumber: _angular_wavenumber,
})

_phase = phase_law.law.rhs.subs({
    phase_law.angular_wavenumber: _angular_wavenumber,
    phase_law.position: _position,
    phase_law.angular_frequency: _angular_frequency,
    phase_law.time: _time,
})

_solution_expr = law.rhs.subs(wave_phase, _phase)

_eqn = wave_eqn.definition.subs({
    wave_eqn.position: _position,
    wave_eqn.time: _time,
    wave_eqn.phase_speed: _phase_velocity,
})

_eqn_lhs_subs = _eqn.lhs.subs(wave_eqn.displacement(_position, _time), _solution_expr).doit()

_eqn_rhs_subs = _eqn.rhs.subs(wave_eqn.displacement(_position, _time), _solution_expr).doit()

assert expr_equals(_eqn_lhs_subs, _eqn_rhs_subs)


@validate_input(wave_phase_=wave_phase)
@validate_output_same("amplitude_")
def calculate_displacement(
    amplitude_: Quantity,
    wave_phase_: Quantity | float,
) -> Quantity:
    # Solution to the wave equation in the form `A * cos(phi)`

    wave_phase_ = scale_factor(wave_phase_)
    result = law.rhs.subs(
        solution(wave_phase),
        amplitude_ * cos(wave_phase_),
    )
    return Quantity(result)
