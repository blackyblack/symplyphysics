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
## Heat capacity is specific heat capacity (heat capacity per unit mass) times mass.

# Definition: C = c * m
## C - heat capacity
## c - specific heat capacity
## m - mass

# Note
## Heat capacity is an extensive quantity, i.e. its scales linearly with the size of the system,
## whereas specific heat capacity is intensive and does not depend on the size of the system or the amount
## of substance in consideration.

heat_capacity = Symbol("heat_capacity", units.energy / units.temperature)
specific_heat_capacity = Symbol("specific_heat_capacity", units.energy / (units.temperature * units.mass))
mass = symbols.basic.mass

definition = Eq(heat_capacity, specific_heat_capacity * mass)


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    specific_heat_capacity_=specific_heat_capacity,
    mass_=mass,
)
@validate_output(heat_capacity)
def calculate_heat_capacity(
    specific_heat_capacity_: Quantity,
    mass_: Quantity,
) -> Quantity:
    result = definition.rhs.subs({
        specific_heat_capacity: specific_heat_capacity_,
        mass: mass_,
    })
    return Quantity(result)
