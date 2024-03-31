from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

# Description
## The second law of thermodynamics applied to a closed system and an idealized, reversible or quasistatic process
## states that in such a process of transfer of energy as heat to closed thermodynamic system (B), which
## allows energy but not matter exchanges, from an auxiliary thermodynamic system (A), an infinitesimal
## increment in the entropy of system B is defined to result from an infinitesimal transfer of heat to
## system B divided by the common thermodynamic temperature of both system A and system B.

# Law: dS = δQ / T
## dS - infinitesimal increase in entropy of system of interest (B)
## δQ - infinitesimal amount of heat transfered to system of interest (B)
## T - common thermodynamic temperature of the systems (A and B)

# Notes
## - dS is an exact differential (since entropy S is a state function) whereas δQ is an inexact differential
## - Also applicable to actually possible quasistatic irreversible processes without composition change

infinitesimal_entropy_change = Symbol("infinitesimal_entropy_change", units.energy / units.temperature)
infinitesimal_transfer_of_heat = Symbol("infinitesimal_transfer_of_heat", units.energy)
common_temperature = clone_symbol(symbols.thermodynamics.temperature, "common_temperature")

law = Eq(infinitesimal_entropy_change, infinitesimal_transfer_of_heat / common_temperature)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    infinitesimal_transfer_of_heat_=infinitesimal_transfer_of_heat,
    common_temperature_=common_temperature,
)
@validate_output(infinitesimal_entropy_change)
def calculate_infinitesimal_entropy_change(
    infinitesimal_transfer_of_heat_: Quantity,
    common_temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        infinitesimal_transfer_of_heat: infinitesimal_transfer_of_heat_,
        common_temperature: common_temperature_,
    })
    return Quantity(result)
