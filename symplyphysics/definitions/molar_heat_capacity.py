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
## Heat capacity is molar heat capacity (heat capacity per unit amount of substance) times amount of substance.

# Definition: C = C_m * n
## C - heat capacity
## C_m - molar heat capacity
## n - amount of substance

# Note
## Heat capacity is an extensive quantity, i.e. its scales linearly with the size of the system,
## whereas molar heat capacity is intensive and does not depend on the size of the system or the amount
## of substance in consideration.

heat_capacity = Symbol("heat_capacity", units.energy / units.temperature)
molar_heat_capacity = Symbol("molar_heat_capacity", units.energy / (units.temperature * units.amount_of_substance))
amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)

definition = Eq(heat_capacity, molar_heat_capacity * amount_of_substance)


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    molar_heat_capacity_=molar_heat_capacity,
    amount_of_substance_=amount_of_substance,
)
@validate_output(heat_capacity)
def calculate_heat_capacity(
    molar_heat_capacity_: Quantity,
    amount_of_substance_: Quantity,
) -> Quantity:
    result = definition.rhs.subs({
        molar_heat_capacity: molar_heat_capacity_,
        amount_of_substance: amount_of_substance_,
    })
    return Quantity(result)
