from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Admittance (a.k.a. complex conductance) is ability of dipole to conduct electrical signal.

# Definition: Y = 1 / Z, where
## Y is admittance of dipole,
## Z is its impedance.

dipole_admittance = Symbol("dipole_admittance", units.conductance)
dipole_impedance = Symbol("dipole_impedance", units.impedance)

definition = Eq(dipole_admittance, 1 / dipole_impedance)

definition_units_SI = units.siemens


def print_law() -> str:
    return print_expression(definition)


@validate_input(impedance_=dipole_impedance)
@validate_output(dipole_admittance)
def calculate_admittance(impedance_: Quantity) -> Quantity:
    solved = solve(definition, dipole_admittance, dict=True)[0][dipole_admittance]
    result_expr = solved.subs({dipole_impedance: impedance_})
    return Quantity(result_expr)
