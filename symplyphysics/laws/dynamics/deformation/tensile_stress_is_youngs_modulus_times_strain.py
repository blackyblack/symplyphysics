from sympy import Eq, Abs
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## When an object is under stress or compression, the stress is related to the strain via the
## Young's modulus.

# Law: sigma = E * |dL|/L
## sigma - stress
## E - Young's modulus for the object
## dL - change in length of object
## L - initial length of object

# Note
## |dL|/L is called the tensile/compressive strain of the object

stress = Symbol("stress", units.pressure)
youngs_modulus = Symbol("youngs_modulus", units.pressure)
change_in_length = Symbol("change_in_length", units.length)
initial_length = Symbol("initial_length", units.length)

law = Eq(stress, youngs_modulus * Abs(change_in_length) / initial_length)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    youngs_modulus_=youngs_modulus,
    change_in_length_=change_in_length,
    initial_length_=initial_length,
)
@validate_output(stress)
def calculate_tensile_stress(
    youngs_modulus_: Quantity,
    change_in_length_: Quantity,
    initial_length_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        youngs_modulus: youngs_modulus_,
        change_in_length: change_in_length_,
        initial_length: initial_length_,
    })
    return Quantity(result)
