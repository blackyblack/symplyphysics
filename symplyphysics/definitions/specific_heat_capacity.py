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
## Specific heat capacity of a substance is its heat capacity divided by its mass.

# Definition: c = C / m
## c - specific heat capacity
## C - heat capacity
## m - mass

# Note
## Heat capacity is an extensive quantity, i.e. its scales linearly with the size of the system,
## whereas specific heat capacity is intensive and does not depend on the size of the system or the amount
## of substance in consideration.

specific_heat_capacity = Symbol("specific_heat_capacity", units.energy / (units.temperature * units.mass))
heat_capacity = Symbol("heat_capacity", units.energy / units.temperature)
mass = symbols.basic.mass

definition = Eq(specific_heat_capacity, heat_capacity / mass)


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    heat_capacity_=heat_capacity,
    mass_=mass,
)
@validate_output(specific_heat_capacity)
def calculate_specific_heat_capacity(
    heat_capacity_: Quantity,
    mass_: Quantity,
) -> Quantity:
    result = definition.rhs.subs({
        heat_capacity: heat_capacity_,
        mass: mass_,
    })
    return Quantity(result)
