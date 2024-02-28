from sympy import Eq, pi, sin
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
from symplyphysics.core.symbols.quantities import scale_factor

# Description
## When an oscillating external force is driving the oscillations of an oscillator,
## the amplitude of oscillations is greatest when the angular frequency of the driving
## force is equal to the natural frequency of the oscillator. This condition is called
## resonance.

# Law: q(t) = f_m * (w0 * t * sin(w0 * t + phi) - sin(phi) * sin(w0*t)) / (2 * m * w0**2)
## q(t) - displacement of resonant oscillations
## t - time
## m - mass of oscillator
## w0 - natural frequency of oscillator
## f_m - amplitude of driving force
## phi - phase lag of driving force

# Conditions
## - Angular frequency of the driving force is equal to the natural frequency of the oscillator.
## - No damping is present in the system.

# Note
## - The expression of the driving force has the form `f_m * cos(w * t + phi)` where `w` is
##   the angular frequency of its oscillations.

resonant_displacement = Function("resonant_displacement", units.length)
oscillator_mass = Symbol("oscillator_mass", units.mass)
natural_angular_frequency = Symbol("natural_angular_frequency", angle_type / units.time)
driving_force_amplitude = Symbol("driving_force_amplitude", units.force)
driving_phase_lag = Symbol("driving_phase_lag", angle_type)
time = Symbol("time", units.time)

law = Eq(
    resonant_displacement(time),
    driving_force_amplitude
    * (natural_angular_frequency * time * sin(natural_angular_frequency * time + driving_phase_lag)
       - sin(driving_phase_lag) * sin(natural_angular_frequency * time))
    / (2 * oscillator_mass * natural_angular_frequency**2)
)

# TODO: derive law from driven oscillations equation


def print_law() -> str:
    return print_expression(law)


@validate_input(
    oscillator_mass_=oscillator_mass,
    natural_angular_frequency_=natural_angular_frequency,
    driving_force_amplitude_=driving_force_amplitude,
    driving_phase_lag_=driving_phase_lag,
    time_=time,
)
@validate_output(resonant_displacement)
def calculate_resonant_displacement(
    oscillator_mass_: Quantity,
    natural_angular_frequency_: Quantity,
    driving_force_amplitude_: Quantity,
    driving_phase_lag_: Quantity | float,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        oscillator_mass: oscillator_mass_,
        natural_angular_frequency: natural_angular_frequency_,
        driving_force_amplitude: driving_force_amplitude_,
        driving_phase_lag: scale_factor(driving_phase_lag_),
        time: time_,
    })
    return Quantity(result)
