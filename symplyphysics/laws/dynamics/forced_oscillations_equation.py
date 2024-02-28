from sympy import Derivative, Eq, cos, dsolve
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
## Forced, or driven, oscillations are a type of oscillations in the precence of an external driving
## force acting on the oscillating system. In the case of an oscillating external force, two angular 
## frequencies are associated with such a system: (1) the natural angular frequency of the system, 
## which is the angular frequency the system would oscillate with if no external force were present, 
## and (2) the angular frequency of the external force driving the oscillations. Such systems can
## undergo resonance if the angular frequency of the driving force is close to the natural angular
## frequency of the oscillator.

# Law: d**2(q(t))/dt**2 + omega0**2*q(t) = (f_m / m) * cos(omega * t + phi)
## q(t) - displacement of the oscillating body
## t - time
## omega0 - natural frequency of the oscillator
## m - mass of the oscillating body
## f_m - amplitude of the driving force
## omega - angular frequency of the driving force
## phi - phase lag of the driving force

displacement = Function("displacement", units.length)
oscillator_mass = Symbol("oscillator_mass", units.mass)
natural_angular_frequency = Symbol("natural_angular_frequency", angle_type / units.time)
driving_force_amplitude = Symbol("driving_force_amplitude", units.force)
driving_angular_frequency = Symbol("driving_angular_frequency", angle_type / units.time)
driving_phase_lag = Symbol("driving_phase_lag", angle_type)
time = Symbol("time", units.time)

law = Eq(
    Derivative(displacement(time), time, 2) + natural_angular_frequency**2 * displacement(time),
    (driving_force_amplitude / oscillator_mass) * cos(driving_angular_frequency * time + driving_phase_lag),
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    initial_position_=displacement,
    initial_velocity_=units.velocity,
    oscillator_mass_=oscillator_mass,
    natural_angular_frequency_=natural_angular_frequency,
    driving_force_amplitude_=driving_force_amplitude,
    driving_angular_frequency_=driving_angular_frequency,
    driving_phase_lag_=driving_phase_lag,
    time_=time,
)
@validate_output(displacement)
def calculate_displacement(
    initial_position_: Quantity,
    initial_velocity_: Quantity,
    oscillator_mass_: Quantity,
    natural_angular_frequency_: Quantity,
    driving_force_amplitude_: Quantity,
    driving_angular_frequency_: Quantity,
    driving_phase_lag_: Quantity,
    time_: Quantity,
) -> Quantity:
    initial_conditions = {
        displacement(0): initial_position_,
        displacement(time).diff(time).subs(time, 0): initial_velocity_,
    }
    eqn = (law.rhs - law.lhs).subs({
        natural_angular_frequency: natural_angular_frequency_,
        driving_angular_frequency: driving_angular_frequency_,
    })
    dsolved = dsolve(eqn, displacement(time), ics=initial_conditions).rhs
    result = dsolved.subs({
        oscillator_mass: oscillator_mass_,
        driving_force_amplitude: driving_force_amplitude_,
        driving_phase_lag: driving_phase_lag_,
        time: time_,
    })
    return Quantity(result)
