"""
Phase shift between two points
==============================

The phase shift of wave oscillations between two points is proportional to the ratio
of the distance between the points and the wavelength of the wave.

**Conditions:**

#. The wave is moving in a one-dimensional environment.

**Links:**

#. `BYJU's <https://byjus.com/physics/relation-between-phase-difference-and-path-difference/>`__.
"""

from sympy import Eq, solve, pi
from symplyphysics import (Quantity, validate_input, validate_output, convert_to_float, symbols,
    clone_as_symbol)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.definitions import angular_wavenumber_is_inverse_wavelength as _wavenumber_def
from symplyphysics.laws.waves import phase_of_traveling_wave as _phase_law

phase_shift = symbols.phase_shift
"""
:symbols:`phase_shift` between the two points.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between the two points.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the wave.
"""

law = Eq(phase_shift, 2 * pi * distance / wavelength)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law for a one-dimensional wave.

_first_position = clone_as_symbol(symbols.position, subscript="1")
_second_position = clone_as_symbol(symbols.position, subscript="2")

# Phase is measured in both points at the same time
_first_phase = _phase_law.law.rhs.subs(_phase_law.position, _first_position)
_second_phase = _phase_law.law.rhs.subs(_phase_law.position, _second_position)

# We can also subtract the first point's phase from the second point's phase, but then we should
# also flip the sign in the RHS of the distance equation.
_phase_shift_eqn = Eq(phase_shift, _first_phase - _second_phase)
_distance_eqn = Eq(distance, _first_position - _second_position)

_wavenumber_eqn = _wavenumber_def.definition.subs({
    _wavenumber_def.angular_wavenumber: _phase_law.angular_wavenumber,
    _wavenumber_def.wavelength: wavelength,
})

_phase_shift_derived = solve(
    (_phase_shift_eqn, _distance_eqn, _wavenumber_eqn),
    (phase_shift, _first_position, _phase_law.angular_wavenumber),
    dict=True,
)[0][phase_shift]

assert expr_equals(_phase_shift_derived, law.rhs)


@validate_input(distance_between_points_=distance, wavelength_=wavelength)
@validate_output(phase_shift)
def calculate_phase_difference(distance_between_points_: Quantity, wavelength_: Quantity) -> float:
    result_expr = solve(law, phase_shift, dict=True)[0][phase_shift]
    result_expr = result_expr.subs({
        distance: distance_between_points_,
        wavelength: wavelength_
    }).doit()
    return convert_to_float(result_expr)


# UNIQUE_LAW_ID: 385
