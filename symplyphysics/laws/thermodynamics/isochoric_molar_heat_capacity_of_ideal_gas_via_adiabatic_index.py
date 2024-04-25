from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## For ideal gases the isobaric and isochoric heat capacities can be calculated via the molar gas
## constant and the adiabatic index of the gas.

# Law: C_V = R / (gamma - 1)
## C_V - isochoric molar heat capacity
## R - molar gas constant
## gamma - adiabatic index, or [heat capacity ratio](../../definitions/heat_capacity_ratio.py)

isochoric_molar_heat_capacity = Symbol(
    "isochoric_molar_heat_capacity",
    units.energy / (units.temperature * units.amount_of_substance)
)
adiabatic_index = Symbol("adiabatic_index", dimensionless)

law = Eq(isochoric_molar_heat_capacity, units.molar_gas_constant / (adiabatic_index - 1))


@validate_input(adiabatic_index_=adiabatic_index)
@validate_output(isochoric_molar_heat_capacity)
def calculate_isochoric_molar_heat_capacity(adiabatic_index_: float) -> Quantity:
    result = law.rhs.subs(adiabatic_index, adiabatic_index_)
    return Quantity(result)
