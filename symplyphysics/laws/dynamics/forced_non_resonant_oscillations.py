from sympy import Eq, cos
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    angle_type,
    validate_input,
    validate_output,
    print_expression,
)
from symplyphysics.core.symbols.quantities import scale_factor

# Description
## Forced, or driven, oscillations are a type of oscillations in the precence of an external driving
## force acting on the oscillating system. In the case of an oscillating external force, two angular 
## frequencies are associated with such a system: (1) the natural angular frequency of the system, 
## which is the angular frequency the system would oscillate with if no external force were present, 
## and (2) the angular frequency of the external force driving the oscillations.

# Law: q(t) = (f / m) * cos(omega*t + phi) / (omega0**2 - omega**2)
## q(t) - particular solution of the forced oscillations equation that accounts 
##        for the oscillator's response to the driving force
## t - time
## m - mass of oscillating body
## omega0 - natural angular frequency of the oscillator
## f - amplitudue of the external driving force
## omega - angular frequency of the external driving force
## phi - phase lag of the oscillations of the external force

# Conditions
## - Angular frequency of the external force is strictly not equal to the natural angular 
##   frequency of the oscillator.
## - No damping is present in the system.

# Notes
## - The external driving force has the form of f(t) = f_m * cos(omega*t + phi).
## - The complete expression of the displacement function can be found as the sum of
##   the solution of simple harmonic motion equation and the particular solution presented here.

driven_displacement = Function("driven_displacement", units.length)
oscillator_mass = Symbol("oscillator_mass", units.mass)
natural_angular_frequency = Symbol("natural_angular_frequency", angle_type / units.time)
driving_force_amplitude = Symbol("driving_force_amplitude", units.force)
driving_angular_frequency = Symbol("driving_angular_frequency", angle_type / units.time)
driving_phase_lag = Symbol("driving_phase_lag", angle_type)
time = Symbol("time", units.time)

law = Eq(
    driven_displacement(time),
    (driving_force_amplitude / oscillator_mass)
    * cos(driving_angular_frequency * time + driving_phase_lag)
    / (natural_angular_frequency**2 - driving_angular_frequency**2)
)

# TODO: prove this is a solution of the driven oscillations equation


def print_law() -> str:
    return print_expression(law)


@validate_input(
    oscillator_mass_=oscillator_mass,
    natural_angular_frequency_=natural_angular_frequency,
    driving_force_amplitude_=driving_force_amplitude,
    driving_angular_frequency_=driving_angular_frequency,
    driving_phase_lag_=driving_phase_lag,
    time_=time,
)
@validate_output(driven_displacement)
def calculate_driven_displacement(
    oscillator_mass_: Quantity,
    natural_angular_frequency_: Quantity,
    driving_force_amplitude_: Quantity,
    driving_angular_frequency_: Quantity,
    driving_phase_lag_: Quantity | float,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        oscillator_mass: oscillator_mass_,
        natural_angular_frequency: natural_angular_frequency_,
        driving_force_amplitude: driving_force_amplitude_,
        driving_angular_frequency: driving_angular_frequency_,
        driving_phase_lag: scale_factor(driving_phase_lag_),
        time: time_,
    })
    return Quantity(result)
