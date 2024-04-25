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

# Law: C_p = R * gamma / (gamma - 1)
## C_p - isobaric molar heat capacity
## R - molar gas constant
## gamma - adiabatic index, or [heat capacity ratio](../../definitions/heat_capacity_ratio.py)

isobaric_molar_heat_capacity = Symbol(
    "isobaric_molar_heat_capacity",
    units.energy / (units.temperature * units.amount_of_substance)
)
adiabatic_index = Symbol("adiabatic_index", dimensionless)

law = Eq(
    isobaric_molar_heat_capacity,
    units.molar_gas_constant * adiabatic_index / (adiabatic_index - 1)
)


@validate_input(adiabatic_index_=adiabatic_index)
@validate_output(isobaric_molar_heat_capacity)
def calculate_isobaric_molar_heat_capacity(adiabatic_index_: float) -> Quantity:
    result = law.rhs.subs(adiabatic_index, adiabatic_index_)
    return Quantity(result)
