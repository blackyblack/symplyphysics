from sympy import Eq, Integral
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    symbols,
)

# Description
## If the equation of state is known, the internal energy of a substance can be found
## as a function of volume at constant temperature.

# Law: u = Integral(c_V(T), T) - a / v
## u - molar internal energy of van der Waals fluid at constant temperature
## c_V - isochoric molar heat capacity as a function of temperature
## T - absolute temperature
## a - parameter of van der Waals equation of state corresponding to magnitude of intermolecular forces
## v - molar volume

# Conditions
## - The fluid is homogeneous and in a single phase state.

molar_internal_energy = Symbol("molar_internal_energy", units.energy / units.amount_of_substance)
isochoric_molar_heat_capacity = Function(
    "isochoric_molar_heat_capacity",
    units.energy / (units.temperature * units.amount_of_substance),
)
temperature = symbols.thermodynamics.temperature
bonding_forces_parameter = Symbol(
    "bonding_forces_parameter",
    units.pressure * (units.volume / units.amount_of_substance)**2,
)
molar_volume = Symbol("molar_volume", units.volume / units.amount_of_substance)

law = Eq(
    molar_internal_energy,
    Integral(isochoric_molar_heat_capacity(temperature), temperature) -
    bonding_forces_parameter / molar_volume,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    isochoric_molar_heat_capacity_=isochoric_molar_heat_capacity,
    temperature_=temperature,
    bonding_forces_parameter_=bonding_forces_parameter,
    molar_volume_=molar_volume,
)
@validate_output(molar_internal_energy)
def calculate_internal_energy(
    isochoric_molar_heat_capacity_: Quantity,
    temperature_: Quantity,
    bonding_forces_parameter_: Quantity,
    molar_volume_: Quantity,
) -> Quantity:
    # Note that internal energy is only known up to a constant term
    # Isochoric heat capacity is assumed to be a constant independent of temperature

    isochoric_molar_heat_capacity_function = isochoric_molar_heat_capacity_
    result = law.rhs.subs(isochoric_molar_heat_capacity(temperature),
        isochoric_molar_heat_capacity_function).doit().subs({
        temperature: temperature_,
        bonding_forces_parameter: bonding_forces_parameter_,
        molar_volume: molar_volume_,
        })
    return Quantity(result)
