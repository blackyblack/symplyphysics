from sympy import Eq, pi
from symplyphysics import (
    units,
    angle_type,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## ...

# Law: q(t) = (f_m / m) * t * cos(w0 * t + phi - pi/2) / (2 * w0)
## q(t) - displacement of resonant oscillations
## t - time
## m - mass of oscillator
## w0 - natural frequency of oscillator
## f_m - amplitude of driving force
## phi - phase lag of driving force

resonant_displacement = Function("resonant_displacement", units.length)
oscillator_mass = Symbol("oscillator_mass", units.mass)
natural_angular_frequency = Symbol("natural_angular_frequency", angle_type / units.time)
driving_force_amplitude = Symbol("driving_force_amplitude", units.force)
driving_phase_lag = Symbol("driving_phase_lag", angle_type)
time = Symbol("time")

law = Eq(
    resonant_displacement(time),
    (driving_force_amplitude / oscillator_mass)
    * cos(natural_angular_frequency * time + driving_phase_lag - pi/2)
    / (2 * natural_angular_frequency)
)
