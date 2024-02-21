from sympy import Eq, exp, dsolve, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import damped_harmonic_oscillator_equation as damped_eqn

# Description
## In the presence of a damping force in the oscillating system, the system's behaviour
## depends on the value of the damping ratio. When it is less than 1, the system oscillates
## with a slightly different frequency than the undamped case with the amplitude decreasing
## to zero. This behaviour is also known as underdamping.

# Law: x(t) = x0 * exp(-lambda*t) * cos(w1*t)

# underdamped_magnitude = symbols("underdamped_magnitude", positive=True)
# underdamped_angular_frequency = undamped_angular_frequency * sqrt(1 - damping_ratio**2)
# underdamped_exponential_decay = undamped_angular_frequency * damping_ratio
# underdamped_phase_lag = symbols("underdamped_phase_lag", real=True)

# # The solition in the underdamped case take the form of
# ## x_underdamped(t) = x0 * exp(-lambda*t) * cos(w1*t + phi), where
# ## - x0 is the amplitude of underdamped oscillations,
# ## - lambda is the exponential decay,
# ## - phi is the phase lag

# underdamped_displacement = (
#     underdamped_magnitude
#     * exp(-1 * underdamped_exponential_decay * time)
#     * cos(underdamped_angular_frequency * time + underdamped_phase_lag)
# )

# _underdamped_subs = (
#     definition
#     .subs(displacement(time), underdamped_displacement)
#     .doit()
#     .simplify()
# )
# assert expr_equals(_underdamped_subs, 0)