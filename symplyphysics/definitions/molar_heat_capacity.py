from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Molar heat capacity of a substance is its heat capacity divided by its amount in moles.

# Definition: C_m = C / n
## C_m - molar heat capacity
## C - heat capacity
## n - amount of substance

# Note
## Heat capacity is an extensive quantity, i.e. its scales linearly with the size of the system,
## whereas molar heat capacity is intensive and does not depend on the size of the system or the amount
## of substance in consideration.

molar_heat_capacity = Symbol("molar_heat_capacity", units.energy / (units.temperature * units.amount_of_substance))
heat_capacity = Symbol("heat_capacity", units.energy / units.temperature)
amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)

definition = Eq(molar_heat_capacity, heat_capacity / amount_of_substance)


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    heat_capacity_=heat_capacity,
    amount_of_substance_=amount_of_substance,
)
@validate_output(molar_heat_capacity)
def calculate_molar_heat_capacity(
    heat_capacity_: Quantity,
    amount_of_substance_: Quantity,
) -> Quantity:
    result = definition.rhs.subs({
        heat_capacity: heat_capacity_,
        amount_of_substance: amount_of_substance_,
    })
    return Quantity(result)
