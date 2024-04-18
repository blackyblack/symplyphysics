from sympy import Eq, solve
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The Gibbs—Duhem relation is a relationship among the intensive parameters of the system.
## Subsequently, For a system with `I` components, there is `I + 1` independent parameters,
## or degrees of freedom.

# Law: S * dT - V * dp + N * d(mu) = 0
## S - entropy
## T - absolute temperature
## V - volume
## p - pressure
## N - particle count
## mu - chemical potential
## Notation: d(x) - exact differential of x

entropy = Symbol("entropy", units.energy / units.temperature)
temperature_change = Symbol("temperature_change", units.temperature)
volume = Symbol("volume", units.volume)
pressure_change = Symbol("pressure_change", units.pressure)
particle_count = Symbol("particle_count", dimensionless)
chemical_potential_change = Symbol("chemical_potential_change", units.energy)

law = Eq(
    entropy * temperature_change - volume * pressure_change + particle_count * chemical_potential_change,
    0,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    entropy_=entropy,
    temperature_change_=temperature_change,
    volume_=volume,
    pressure_change_=pressure_change,
    particle_count_=particle_count,
)
@validate_output(chemical_potential_change)
def calculate_chemical_potential_change(
    entropy_: Quantity,
    temperature_change_: Quantity,
    volume_: Quantity,
    pressure_change_: Quantity,
    particle_count_: int,
) -> Quantity:
    result = solve(law, chemical_potential_change)[0].subs({
        entropy: entropy_,
        temperature_change: temperature_change_,
        volume: volume_,
        pressure_change: pressure_change_,
        particle_count: particle_count_,
    })
    return Quantity(result)
