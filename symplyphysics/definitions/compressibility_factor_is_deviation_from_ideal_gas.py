from sympy import Eq, S
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
    convert_to,
)

# Description
## The compressibility factor, also known as the compression factor or the gas deviation factor,
## describes the deviation of a real gas from ideal gas behaviour. In general, the deviation from
## ideal gas behaviour becomes more prominent the closer the gas is to a phase change, the lower
## the temperature or the larger the pressure.

# Definition: Z = (p * V) / (n * R * T)
## Z - compressibility factor
## p - pressure
## V - volume
## n - amount of substance
## R - molar gas constant
## T - temperature

# Note
## - Can be equivalently defined as the ratio of the molar volume of the real gas (V / n) to the
##   molar volume of the corresponding ideal gas (R * T / p) at the same temperature and pressure.
## - Z = 1 is the case of ideal gas behaviour.
## - At high pressures molecules collide more often leading to an increase of repulsive forces between
##   molecules, making the molar volume of the real gas (V / n) greater than that of ideal gas (R * T / p),
##   in other words the particles have a larger extended volume, leading to Z > 1.
## - At lower pressures, molecules are free to move and attractive forces dominate, leading to Z < 1.

compressibility_factor = Symbol("compressibility_factor", dimensionless)
pressure = Symbol("pressure", units.pressure)
volume = Symbol("volume", units.volume)
amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)
temperature = symbols.thermodynamics.temperature

definition = Eq(
    compressibility_factor,
    (pressure * volume) / (amount_of_substance * units.molar_gas_constant * temperature)
)


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    pressure_=pressure,
    volume_=volume,
    amount_of_substance_=amount_of_substance,
    temperature_=temperature,
)
@validate_output(compressibility_factor)
def calculate_compressibility_factor(
    pressure_: Quantity,
    volume_: Quantity,
    amount_of_substance_: Quantity,
    temperature_: Quantity,
) -> float:
    result = definition.rhs.subs({
        pressure: pressure_,
        volume: volume_,
        amount_of_substance: amount_of_substance_,
        temperature: temperature_,
    })
    return float(convert_to(Quantity(result), S.One))
