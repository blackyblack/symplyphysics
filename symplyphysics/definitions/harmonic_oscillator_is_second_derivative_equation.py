from sympy import (Derivative, Eq, cos, solve, symbols, Function as SymFunction)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals

# Description
## In classical mechanics, a harmonic oscillator is a system that, when displaced from its equilibrium position, experiences a restoring force F proportional to the displacement x.
## If F is the only force acting on the system, the system is called a simple harmonic oscillator.
## Displacement is not only limited to physical motion, but should be interpreted in general terms. Harmonic oscillator can represent mechanical systems that include pendulums
## (with small angles of displacement), masses connected to springs, and acoustical systems. Other analogous systems include electrical harmonic oscillators such as RLC circuits.

# Definition: d^2x(t)/dt^2 = -w**2 * x(t)
# Where:
## x(t) is function of dispalacement of the oscillator over time,
## t is time,
## w is angular frequency,
##   see [angular frequency](../laws/kinematic/period_from_angular_frequency.py) implementation,
## d^2 is second derivative.

# Conditions:
## - No frictional force (damping).

# Can have various dimensions - do not specify it
displacement_function = symbols("displacement", cls=SymFunction)
angular_frequency = Symbol("angular_frequency", units.frequency)
time = Symbol("time", units.time)

definition = Eq(Derivative(displacement_function(time), (time, 2)),
    -angular_frequency**2 * displacement_function(time))

# Confirm that cosine function is a solution to this equation

## Expected solution for x"(t) = -w**2 * x(t) is:
## A * e^(i * w * t) + B * e^(-i * w * t)
## This form is known to be represented as cosine function:
## A * cos(w * t + phi), or sine function:
## A * sin(w * t + phi + pi / 2)

## dsolve() gives us solution in exponential form - we are looking for the solution as trigonometric function

amplitude = symbols("amplitude")
initial_phase = symbols("initial_phase")
displacement_function_eq = Eq(displacement_function(time),
    amplitude * cos(angular_frequency * time + initial_phase))
dsolved = definition.subs(displacement_function(time), displacement_function_eq.rhs)
assert expr_equals(dsolved.lhs, dsolved.rhs)

## There are many solutions for harmonic_oscillation_eq. Add condition, that at initial point of time (time = 0)
## there is max displacement (displacement(time) = amplitude).
## Let's prove that initial phase of cosine function (displacement_function_eq) should be zero.

initial_condition = Eq(displacement_function(0), amplitude)
displacement_function_at_zero_time_eq = displacement_function_eq.subs(time, 0)
## Initial phase solutions have period of 2*pi. Take first solution.
initial_phase_solved = solve([displacement_function_at_zero_time_eq, initial_condition],
    (amplitude, initial_phase),
    dict=True)[0][initial_phase]
assert expr_equals(initial_phase_solved, 0)


def print_law() -> str:
    return print_expression(definition)


@validate_input(amplitude_=units.length, angular_frequency_=angular_frequency, time_=time)
@validate_output(units.length)
def calculate_displacement(amplitude_: Quantity, angular_frequency_: Quantity,
    time_: Quantity) -> Quantity:
    result_expr = displacement_function_eq.subs({
        amplitude: amplitude_,
        initial_phase: 0,
        angular_frequency: angular_frequency_,
        time: time_
    }).rhs
    return expr_to_quantity(result_expr)
