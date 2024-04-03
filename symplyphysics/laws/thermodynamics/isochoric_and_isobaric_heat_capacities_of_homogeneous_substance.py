from sympy import Eq, solve
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
## Mayer's relation is the relation between heat capacity at constant pressure and heat
## capacity at constant volume.

# Law: C_p - C_V = V * T * (alpha_V)**2 / beta_T
## C_p - heat capacity at constant pressure
## C_V - heat capacity at constant volume
## V - volume
## T - temperature
## alpha_V - [thermal expansion coefficient](../../definitions/volumetric_coefficient_of_thermal_expansion.py)
## beta_T - isothermal compressibility, see [this](../../definitions/thermodynamic_compressibility.py) 
##          with the derivative taken at constant temperature

# Note
## - Applicable to any homogeneous substances, not just ideal gases.

isobaric_heat_capacity = Symbol("isobaric_heat_capacity", units.energy / units.temperature)
isochoric_heat_capacity = Symbol("isochoric_heat_capacity", units.energy / units.temperature)
volume = Symbol("volume", units.volume)
temperature = symbols.thermodynamics.temperature
thermal_expansion_coefficient = Symbol("thermal_expansion_coefficient", 1 / units.temperature)
isothermal_compressibility = Symbol("isothermal_compressibility", 1 / units.pressure)

law = Eq(
    isobaric_heat_capacity - isochoric_heat_capacity,
    volume * temperature * thermal_expansion_coefficient**2 / isothermal_compressibility,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    isobaric_heat_capacity_=isobaric_heat_capacity,
    volume_=volume,
    temperature_=temperature,
    thermal_expansion_coefficient_=thermal_expansion_coefficient,
    isothermal_compressibility_=isothermal_compressibility,
)
@validate_output(isochoric_heat_capacity)
def calculate_isochoric_heat_capacity(
    isobaric_heat_capacity_: Quantity,
    volume_: Quantity,
    temperature_: Quantity,
    thermal_expansion_coefficient_: Quantity,
    isothermal_compressibility_: Quantity,
) -> Quantity:
    expr = solve(law, isochoric_heat_capacity)[0]
    result = expr.subs({
        isobaric_heat_capacity: isobaric_heat_capacity_,
        volume: volume_,
        temperature: temperature_,
        thermal_expansion_coefficient: thermal_expansion_coefficient_,
        isothermal_compressibility: isothermal_compressibility_,
    })
    return Quantity(result)
