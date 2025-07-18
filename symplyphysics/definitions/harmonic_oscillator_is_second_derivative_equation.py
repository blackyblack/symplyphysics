"""
Harmonic oscillator is a second order derivative equation
=========================================================

In classical mechanics a *simple harmonic oscillator* is a system that, after a small displacement from equilibrium, experiences a restoring force :math:`F` proportional to that displacement. Displacement is not only limited to physical motion, but should be interpreted in general terms. Examples include small-angle pendulums, mass–spring systems, acoustic resonators, and electrical RLC circuits.
If :math:`F` is the only force acting on the system, the system is called a simple harmonic oscillator.

**Conditions:**

#. There is no damping (i.e. friction) in the system.
#. The system experiences a single restoring force :math:`F` (for mechanical oscillators).

**Links:**

#. `Wikipedia – Simple harmonic oscillator <https://en.wikipedia.org/wiki/Harmonic_oscillator#Simple_harmonic_oscillator>`__
"""

from sympy import Derivative, Eq, cos, solve, symbols as sympy_symbols
from symplyphysics import (
    clone_as_function,
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals

time = symbols.time
"""
:symbols:`time`.
"""

# Can have various dimensions - do not specify it
displacement = clone_as_function(symbols.any_quantity, [time],
    display_symbol="x",
    display_latex="x")
"""
Displacement of oscillator from equilibrium as a function of time. See :symbols:`any_quantity`.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the oscillator.
"""

definition = Eq(Derivative(displacement(time), (time, 2)),
    -angular_frequency**2 * displacement(time))
"""
:laws:symbol::

:laws:latex::
"""

# Confirm that cosine function is a solution to this equation

## Expected solution for x"(t) = -w**2 * x(t) is:
## A * e^(i * w * t) + B * e^(-i * w * t)
## This form is known to be represented as cosine function:
## A * cos(w * t + phi), or sine function:
## A * sin(w * t + phi + pi / 2)

## dsolve() gives us solution in exponential form - we are looking for the solution as trigonometric function

_amplitude = sympy_symbols("amplitude")
_initial_phase = sympy_symbols("initial_phase")
_displacement_function_eq = Eq(displacement(time),
    _amplitude * cos(angular_frequency * time + _initial_phase))
_dsolved = definition.subs(displacement(time), _displacement_function_eq.rhs)
assert expr_equals(_dsolved.lhs, _dsolved.rhs)

## There are many solutions for harmonic_oscillation_eq. Add condition, that at initial point of time (time = 0)
## there is max displacement (displacement(time) = _amplitude).
## Let's prove that initial phase of cosine function (_displacement_function_eq) should be zero.

_initial_condition = Eq(displacement(0), _amplitude)
_displacement_function_at_zero_time_eq = _displacement_function_eq.subs(time, 0)
## Initial phase solutions have period of 2*pi. Take first solution.
_initial_phase_solved = solve([_displacement_function_at_zero_time_eq, _initial_condition],
    (_amplitude, _initial_phase),
    dict=True)[0][_initial_phase]
assert expr_equals(_initial_phase_solved, 0)


@validate_input(amplitude_=units.length, angular_frequency_=angular_frequency, time_=time)
@validate_output(units.length)
def calculate_displacement(amplitude_: Quantity, angular_frequency_: Quantity,
    time_: Quantity) -> Quantity:
    result_expr = _displacement_function_eq.subs({
        _amplitude: amplitude_,
        _initial_phase: 0,
        angular_frequency: angular_frequency_,
        time: time_
    }).rhs
    return Quantity(result_expr)
