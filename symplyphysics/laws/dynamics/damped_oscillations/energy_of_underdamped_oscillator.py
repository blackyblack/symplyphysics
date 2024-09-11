from sympy import Eq, exp
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)

# Description
## In the presence of a damping force the oscillating system is no longer closed and its energy
## dissipates to the environment. The total energy of the oscillator becomes converted into thermal
## energy. For small values of the damping ratio, the following equation approximately describes
## the total mechanical energy of the underdamped oscillator:

# Law: E ~= 1/2 * m * w0**2 * A**2 * exp(-2*lambda*t)
## E - total mechanical energy of oscillator
## m - mass of oscillating object
## A - maximum amplitude of oscillations
## w0 - undamped (natural) angular frequency of oscillations
## lambda - [exponential decay constant](../../kinematics/damped_oscillations/damping_ratio_from_decay_constant_and_undamped_frequency.py)
## t - time

# Conditions
## - Damping ratio (zeta) is small, i.e. zeta << 1

oscillator_energy = Symbol("oscillator_energy", units.energy)
oscillator_mass = clone_symbol(symbols.mass)
maximum_amplitude = Symbol("maximum_amplitude", units.length)
undamped_angular_frequency = Symbol("undamped_angular_frequency", angle_type / units.time)
exponential_decay_constant = Symbol("exponential_decay_constant", angle_type / units.time)
time = Symbol("time", units.time)

law = Eq(
    oscillator_energy, oscillator_mass * undamped_angular_frequency**2 * maximum_amplitude**2 *
    exp(-2 * exponential_decay_constant * time) / 2)

# TODO Derive from [underdamped oscillations](../../kinematics/damped_oscillations/underdamping.py)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    oscillator_mass_=oscillator_mass,
    maximum_amplitude_=maximum_amplitude,
    undamped_angular_frequency_=undamped_angular_frequency,
    exponential_decay_constant_=exponential_decay_constant,
    time_=time,
)
@validate_output(oscillator_energy)
def calculate_oscillator_energy(
    oscillator_mass_: Quantity,
    maximum_amplitude_: Quantity,
    undamped_angular_frequency_: Quantity,
    exponential_decay_constant_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        oscillator_mass: oscillator_mass_,
        maximum_amplitude: maximum_amplitude_,
        undamped_angular_frequency: undamped_angular_frequency_,
        exponential_decay_constant: exponential_decay_constant_,
        time: time_,
    })
    return Quantity(result)
