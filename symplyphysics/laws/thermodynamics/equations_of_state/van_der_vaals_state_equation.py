from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
# To more accurately describe the behavior of real gases at low temperatures,
# a Van der Waals gas model was created, taking into account the forces of intermolecular interaction.
# In this model, internal energy becomes a function not only of temperature, but also of volume.
#
# The Van der Waals equation is one of the well-known approximate equations of state describing
# the properties of a real gas, having a compact form and taking
# into account the main characteristics of a gas with intermolecular interaction.

# Law: (p + a * (nu / V)^2) * (V - b * nu) = nu * R * T
# Where:
## p - pressure
## T - temperature
## V - volume
## nu - amount of substance
## R - universal gas constant
## a - measure of the bonding forces between molecules,
## b is proportional to the volume of 1 mole of molecules.

pressure = Symbol("pressure", units.pressure)
volume = Symbol("volume", units.volume)
temperature = Symbol("temperature", units.temperature)
amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)

bonding_forces_parameter = Symbol("bonding_forces_parameter", units.pressure * units.volume ** 2)
molecules_volume_parameter = Symbol("molecules_volume_parameter", units.volume)


law = Eq(
    (pressure - bonding_forces_parameter * (amount_of_substance / volume) ** 2) * (volume - bonding_forces_parameter * amount_of_substance),
    amount_of_substance * units.molar_gas_constant * temperature
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    volume_=volume,
    temperature_=temperature,
    amount_of_substance_=amount_of_substance,
    bonding_forces_parameter_=bonding_forces_parameter,
    molecules_volume_parameter_=molecules_volume_parameter,
)
@validate_output(pressure)
def calculate_pressure(
    volume_: Quantity,
    temperature_: Quantity,
    amount_of_substance_: Quantity,
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity
) -> Quantity:
    solved = solve(law, pressure, dict=True)[0][pressure]
    result_expr = solved.subs({
        volume: volume_,
        temperature: temperature_,
        amount_of_substance: amount_of_substance_,
        bonding_forces_parameter: bonding_forces_parameter_,
        molecules_volume_parameter: molecules_volume_parameter_
    })
    return Quantity(result_expr)
