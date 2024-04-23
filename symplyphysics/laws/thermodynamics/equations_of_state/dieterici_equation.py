from sympy import Eq, solve, exp
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
)

# Description
## The Dieterici equation is another type of semi-empirical equations approximating real gases
## along with the more well-known van der Waals equation of state.

# Law: p * (Vm - b) = R * T * exp(-a / (R * T * Vm))
## p - pressure
## Vm - molar volume
## T - temperature
## R - molar gas constant
## a, b - parameters specific to each substance

# Notes
## - Like the van der Waals equation of state, the Dieterici equation is also semi-empirical.
## - It approximates moderate pressures of real gases much better than the van der Waals equation,
##   but it is absolutely inapplicable for large pressures.
## - Works only in the limit `b << Vm`, `a << P * Vm**2`
## - Can be converted to the van der Waals equation in the aforementioned limits and the additional
##   limit `a << R * T * V_m`

pressure = Symbol("pressure", units.pressure)
molar_volume = Symbol("molar_volume", units.volume / units.amount_of_substance)
temperature = symbols.thermodynamics.temperature
bonding_forces_parameter = Symbol(
    "bonding_forces_parameter",
    units.pressure * (units.volume / units.amount_of_substance)**2,
)
molecules_volume_parameter = Symbol(
    "molecules_volume_parameter",
    units.volume / units.amount_of_substance,
)

law = Eq(
    pressure * (molar_volume - molecules_volume_parameter),
    units.molar_gas_constant * temperature * exp(-1 * bonding_forces_parameter /
    (units.molar_gas_constant * temperature * molar_volume)))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    molar_volume_=molar_volume,
    temperature_=temperature,
    bonding_forces_parameter_=bonding_forces_parameter,
    molecules_volume_parameter_=molecules_volume_parameter,
)
@validate_output(pressure)
def calculate_pressure(
    molar_volume_: Quantity,
    temperature_: Quantity,
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
) -> Quantity:
    expr = solve(law, pressure)[0]
    result = expr.subs({
        molar_volume: molar_volume_,
        temperature: temperature_,
        bonding_forces_parameter: bonding_forces_parameter_,
        molecules_volume_parameter: molecules_volume_parameter_,
    })
    return Quantity(result)
