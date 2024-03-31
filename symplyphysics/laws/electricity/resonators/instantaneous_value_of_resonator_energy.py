from sympy import Eq, solve, exp
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, angle_type,)

## Description
## A rectangular resonator consists of metal walls and a material filling it.
## The resonator is capable of storing energy. The instantaneous value of the resonator
## energy depends on its quality factor, initial energy value, time and angular frequency.

## Law is: W = W0 * exp(-w0 * t / Q), where
## W - instantaneous value of resonator energy,
## W0 - initial energy of resonator,
## w0 - angular frequency,
## t - time,
## Q - quality factor of resonator.

instantaneous_energy = Symbol("instantaneous_energy", units.energy)

initial_energy = Symbol("initial_energy", units.energy)
time = Symbol("time", units.time)
angular_frequency = Symbol("angular_frequency", angle_type / units.time)
quality_factor = Symbol("quality_factor", dimensionless)

law = Eq(instantaneous_energy, initial_energy * exp(-angular_frequency * time / quality_factor))


def print_law() -> str:
    return print_expression(law)


@validate_input(initial_energy_=initial_energy, time_=time, angular_frequency_=angular_frequency, quality_factor_=quality_factor)
@validate_output(instantaneous_energy)
def calculate_instantaneous_energy(initial_energy_: Quantity, time_: Quantity, angular_frequency_: Quantity,
    quality_factor_: float) -> Quantity:
    result_velocity_expr = solve(law, instantaneous_energy, dict=True)[0][instantaneous_energy]
    result_expr = result_velocity_expr.subs({
        initial_energy: initial_energy_,
        time: time_,
        angular_frequency: angular_frequency_,
        quality_factor: quality_factor_
    })
    return Quantity(result_expr)
