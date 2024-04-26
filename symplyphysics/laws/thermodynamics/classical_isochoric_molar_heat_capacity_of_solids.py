from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_output,
)

# Description
## The Dulong-Petit law states that the classical expression for the molar specific heat capacity 
## of certain chemical elements is constant for temperatures far from the absolute zero.

# Law: C_V = 3 * R
## C_V - isochoric molar heat capacity
## R - molar gas constant

# Conditions
## - The temperature of the system is big enough to omit quantum effects.

isochoric_molar_heat_capacity = Symbol(
    "isochoric_molar_heat_capacity",
    units.energy / (units.temperature * units.amount_of_substance)
)

law = Eq(isochoric_molar_heat_capacity, 3 * units.molar_gas_constant)


@validate_output(isochoric_molar_heat_capacity)
def calculate_isochoric_molar_heat_capacity() -> Quantity:
    return Quantity(law.rhs)
