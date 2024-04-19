from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Mayer's relation is the relation between heat capacity at constant pressure and that at
## constant volume in the case of an ideal gas.

# Law: C_p - C_V = n * R
## C_p - heat capacity at constant pressure
## C_V - heat capacity at constant volume
## n - amount of substance
## R - molar gas constant

isobaric_heat_capacity = Symbol("isobaric_heat_capacity", units.energy / units.temperature)
isochoric_heat_capacity = Symbol("isochoric_heat_capacity", units.energy / units.temperature)
amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)

law = Eq(
    isobaric_heat_capacity - isochoric_heat_capacity,
    amount_of_substance * units.molar_gas_constant,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    isobaric_heat_capacity_=isobaric_heat_capacity,
    amount_of_substance_=amount_of_substance,
)
@validate_output(isochoric_heat_capacity)
def calculate_isochoric_heat_capacity(
    isobaric_heat_capacity_: Quantity,
    amount_of_substance_: Quantity,
) -> Quantity:
    expr = solve(law, isochoric_heat_capacity)[0]
    result = expr.subs({
        isobaric_heat_capacity: isobaric_heat_capacity_,
        amount_of_substance: amount_of_substance_,
    })
    return Quantity(result)
