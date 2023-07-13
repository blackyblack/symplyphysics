from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input, validate_output)

# Description
## Conductivity is ability of medium to conduct electrical current.

# Definition: sigma = 1 / R
# Where:
## sigma is electrical conductivity of object,
## R is resistance

object_conductivity = Symbol("object_conductivity", units.conductance)
object_resistance = Symbol("object_resistance", units.impedance)

definition = Eq(object_conductivity, 1 / object_resistance)

definition_units_SI = units.siemens


def print_law() -> str:
    return print_expression(definition)


@validate_input(resistance_=object_resistance)
@validate_output(object_conductivity)
def calculate_conductivity(resistance_: Quantity) -> Quantity:
    solved = solve(definition, object_conductivity, dict=True)[0][object_conductivity]
    result_expr = solved.subs({object_resistance: resistance_})
    return Quantity(result_expr)
