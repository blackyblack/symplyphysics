from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    dimensionless, validate_input, validate_output)

# Description
## Volume number density is the number of specified objects per unit volume.

# Definition: N = n / V
# Where:
## n is the total number of objects
## V is volume

number_density = Symbol("number_density", 1 / units.volume)
objects = Symbol("objects", dimensionless)
volume = Symbol("volume", units.volume)

definition = Eq(number_density, objects / volume)

definition_units_SI = 1 / units.meter**3


def print_law() -> str:
    return print_expression(definition)


@validate_input(objects_=objects, volume_=volume)
@validate_output(number_density)
def calculate_number_density(objects_: int, volume_: Quantity) -> Quantity:
    solved = solve(definition, number_density, dict=True)[0][number_density]
    result_expr = solved.subs({objects: objects_, volume: volume_})
    return expr_to_quantity(result_expr)
