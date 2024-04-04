from sympy import Eq
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
## The second virial coefficient is a coefficient appearing in the virial equation of state of
## a gas which describes the deviation from ideal gas behaviour. This coefficient describes pair
## interaction between molecules of the substance and can be found if the pair interaction potential
## is known, but it can also be derived from the equation of state as a series expansion with respect
## to inverse molar volume or equivalently to molar density.

# Law: B = b - a / (R * T)
## B - second virial coefficient of the virial expansion
## a - parameter of the van der Waals equation of state representing the magnitude of intermolecular forces
## b - parameter of the van der Waals equation of state representing the effective molecular size
## R - molar gas constant
## T - absolute temperature

second_virial_coefficient = Symbol("second_virial_coefficient", units.volume / units.amount_of_substance)
bonding_forces_parameter = Symbol(
    "bonding_forces_parameter",
    units.pressure * (units.volume / units.amount_of_substance)**2,
)
molecules_volume_parameter = Symbol(
    "molecules_volume_parameter",
    units.volume / units.amount_of_substance,
)
temperature = symbols.thermodynamics.temperature

law = Eq(
    second_virial_coefficient,
    molecules_volume_parameter - bonding_forces_parameter / (units.molar_gas_constant * temperature)
)

# TODO: Derive from the van der Waals equation of state


def print_law() -> str:
    return print_expression(law)


@validate_input(
    bonding_forces_parameter_=bonding_forces_parameter,
    molecules_volume_parameter_=molecules_volume_parameter,
    temperature_=temperature,
)
@validate_output(second_virial_coefficient)
def calculate_second_virial_coefficient(
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
    temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        bonding_forces_parameter: bonding_forces_parameter_,
        molecules_volume_parameter: molecules_volume_parameter_,
        temperature: temperature_,
    })
    return Quantity(result)
