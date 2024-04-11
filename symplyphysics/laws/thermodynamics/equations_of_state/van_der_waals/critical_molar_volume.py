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
## Critical parameters of the van der Waals equation of state are such value of volume, pressure, and
## temperature at which the isotherm has an inflection point whose tangent at that point is zero, i.e.
## the first and second derivatives of pressure with respect to volume at constant temperature are zero.

# Law: V_c = 3 * b
## V_c - critical molar volume
## b - parameter of van der Waals equation of state

critical_molar_volume = Symbol("critical_molar_volume", units.volume / units.amount_of_substance)
molecules_volume_parameter = Symbol(
    "molecules_volume_parameter",
    units.volume / units.amount_of_substance,
)

law = Eq(critical_molar_volume, 3 * molecules_volume_parameter)


def print_law() -> str:
    return print_expression(law)


@validate_input(molecules_volume_parameter_=molecules_volume_parameter)
@validate_output(critical_molar_volume)
def calculate_critical_molar_volume(molecules_volume_parameter_: Quantity) -> Quantity:
    result = law.rhs.subs(molecules_volume_parameter, molecules_volume_parameter_)
    return Quantity(result)
