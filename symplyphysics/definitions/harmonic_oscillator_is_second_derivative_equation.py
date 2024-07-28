"""
Harmonic oscillator is a second order derivative equation
=========================================================

In classical mechanics, a harmonic oscillator is a system that, when displaced from its equilibrium position, experiences a restoring force F proportional to the displacement x.
If F is the only force acting on the system, the system is called a simple harmonic oscillator.
Displacement is not only limited to physical motion, but should be interpreted in general terms. Harmonic oscillator can represent mechanical systems that include pendulums
(with small angles of displacement), masses connected to springs, and acoustical systems. Other analogous systems include electrical harmonic oscillators such as RLC circuits.

**Conditions:**

#. There is no damping (i.e. friction) in the system.
"""

from sympy import (Derivative, Eq, cos, solve, symbols, Function as SymFunction)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals

# Can have various dimensions - do not specify it
displacement_function = symbols("displacement", cls=SymFunction)
"""
Displacement of oscillator from equilibrium as a function of time.

Symbol:
    :code:`x(t)`
"""

angular_frequency = Symbol("angular_frequency", units.frequency)
r"""
Angular frequency of the oscillator.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

definition = Eq(Derivative(displacement_function(time), (time, 2)),
    -angular_frequency**2 * displacement_function(time))
r"""
:code:`Derivative(x(t), (t, 2)) = -1 * w^2 * x(t)`

Latex:
    .. math::
        \frac{d^2 x(t)}{d t^2} = - \omega^2 x(t)
"""

# Confirm that cosine function is a solution to this equation

## Expected solution for x"(t) = -w**2 * x(t) is:
## A * e^(i * w * t) + B * e^(-i * w * t)
## This form is known to be represented as cosine function:
## A * cos(w * t + phi), or sine function:
## A * sin(w * t + phi + pi / 2)

## dsolve() gives us solution in exponential form - we are looking for the solution as trigonometric function

_amplitude = symbols("amplitude")
_initial_phase = symbols("initial_phase")
_displacement_function_eq = Eq(displacement_function(time),
    _amplitude * cos(angular_frequency * time + _initial_phase))
_dsolved = definition.subs(displacement_function(time), _displacement_function_eq.rhs)
assert expr_equals(_dsolved.lhs, _dsolved.rhs)

## There are many solutions for harmonic_oscillation_eq. Add condition, that at initial point of time (time = 0)
## there is max displacement (displacement(time) = _amplitude).
## Let's prove that initial phase of cosine function (_displacement_function_eq) should be zero.

_initial_condition = Eq(displacement_function(0), _amplitude)
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
