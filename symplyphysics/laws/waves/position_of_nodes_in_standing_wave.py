from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Symbol,
    Quantity,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## In a standing wave, the locations of zero amplitude are called nodes. These locations, however,
## are not arbitrary and are integer multiples of half the wavelength of the standing wave.

# Law: x_node = m * (lambda / 2)
## x_node - position of m-th node
## m - integer
## lambda - wavelength

node_position = Symbol("node_position", units.length, real=True)
integer_factor = Symbol("integer_factor", dimensionless, integer=True)
wavelength = Symbol("wavelength", units.length, positive=True)

law = Eq(node_position, integer_factor * wavelength / 2)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    integer_factor_=integer_factor,
    wavelength_=wavelength,
)
@validate_output(node_position)
def calculate_node_position(
    integer_factor_: int, 
    wavelength_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        integer_factor: integer_factor_,
        wavelength: wavelength_,
    })
    return Quantity(result)
