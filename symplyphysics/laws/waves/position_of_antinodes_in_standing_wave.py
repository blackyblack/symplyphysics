from sympy import Eq, S
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
## In a standing wave, the locations of maximum amplitude are called antinodes. These locations are
## not arbitrary, however, and are integer multiples of half the wavelength shifted by a quarter
## of the wavelength:

# Law: x_antinode = (m + 1/2) * (lambda / 2)
## x_antinode - position of m-th antinode
## m - integer
## lambda - wavelength

antinode_position = Symbol("antinode_position", units.length, real=True)
integer_factor = Symbol("integer_factor", dimensionless, integer=True)
wavelength = Symbol("wavelength", units.length, positive=True)

law = Eq(antinode_position, (integer_factor + S(1)/2) * wavelength / 2)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    integer_factor_=integer_factor,
    wavelength_=wavelength,
)
@validate_output(antinode_position)
def calculate_antinode_position(
    integer_factor_: int,
    wavelength_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        integer_factor: integer_factor_,
        wavelength: wavelength_,
    })
    return Quantity(result)
