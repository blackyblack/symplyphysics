from sympy import (I, Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Impedance is the combination of resistance and reactance (both inductive and capacitive) and is
## a complex number, containing both real and imaginary parts. (The real part of impedance is
## resistance, while the imaginary part is reactance.)

# Definition: Z = R + j * X
# Where:
## Z is impedance,
## R is resistance,
## X is reactance,
## j is imaginary number.

impedance = Symbol("impedance", units.impedance)
resistance = Symbol("resistance", units.impedance)
reactance = Symbol("reactance", units.impedance)

definition = Eq(impedance, resistance + I * reactance)

definition_units_SI = units.ohm


def print_law() -> str:
    return print_expression(definition)


@validate_input(resistance_=resistance, reactance_=reactance)
@validate_output(impedance)
def calculate_impedance_magnitude(resistance_: Quantity, reactance_: Quantity) -> Quantity:
    solved = solve(definition, impedance, dict=True)[0][impedance]
    result_expr = solved.subs({resistance: resistance_, reactance: reactance_})
    result_magnitude = abs(result_expr)
    return Quantity(result_magnitude)
