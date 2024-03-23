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
## Coefficients of thermal expansion describe how the size of an object changes with a change in temperature
## at a constant pressure. In isotropic materials, the volumetric coefficient is three times the linear one.

# Law: beta = 3 * alpha
## beta - volumetric coefficient of thermal expansion
## alpha - linear coefficient of thermal expansion

volumetric_expansion_coefficient = Symbol("volumetric_expansion_coefficient", 1 / units.temperature)
linear_expansion_coefficient = Symbol("linear_expansion_coefficient", 1 / units.temperature)

law = Eq(volumetric_expansion_coefficient, 3 * linear_expansion_coefficient)


def print_law() -> str:
    return print_expression(law)


@validate_input(linear_expansion_coefficient_=linear_expansion_coefficient)
@validate_output(volumetric_expansion_coefficient)
def calculate_volumetric_expansion_coefficient(linear_expansion_coefficient_: Quantity) -> Quantity:
    result = law.rhs.subs(linear_expansion_coefficient, linear_expansion_coefficient_)
    return Quantity(result)
