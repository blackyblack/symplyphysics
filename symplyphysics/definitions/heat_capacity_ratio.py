from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    convert_to_float,
)

# Description
## The heat capacity ratio, also known as the adiabatic index, the ratio of specific heats, or
## the insentropic expansion factor, is the ratio of the heat capacity at constant pressure
## to that of constant volume. The heat capacity ratio is used in the description of thermodynamic
## reversible processes; the speed of sound also depends on this factor.

# Law: gamma = C_p / C_V
## gamma - heat capacity ratio
## C_p - heat capacity at constant pressure
## C_V - heat capacity at constant volume

# Note:
## One can also use specific or molar heat capacities in place of the normal heat capacity.

heat_capacity_ratio = Symbol("heat_capacity_ratio", dimensionless)
isobaric_heat_capacity = Symbol("isobaric_heat_capacity", units.energy / units.temperature)
isochoric_heat_capacity = Symbol("isochoric_heat_capacity", units.energy / units.temperature)

definition = Eq(heat_capacity_ratio, isobaric_heat_capacity / isochoric_heat_capacity)


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    isobaric_heat_capacity_=isobaric_heat_capacity,
    isochoric_heat_capacity_=isochoric_heat_capacity,
)
@validate_output(heat_capacity_ratio)
def calculate_heat_capacity_ratio(
    isobaric_heat_capacity_: Quantity,
    isochoric_heat_capacity_: Quantity,
) -> float:
    result = definition.rhs.subs({
        isobaric_heat_capacity: isobaric_heat_capacity_,
        isochoric_heat_capacity: isochoric_heat_capacity_,
    })
    return convert_to_float(result)
