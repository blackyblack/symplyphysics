from sympy import Eq, solve
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)

# Description
## If a standing wave occurs in a string with fixed ends, the value of the displacement function
## must be zero on both ends of the string. To satisfy this boundary condition, only an integer
## number of half-waves can fit into the total length of the string:

# Law: n * (lambda / 2) = L
## n - positive integer, also known as harmonic number of the n-th harmonic
## lambda - wavelength of standing wave
## L - length of the string

integer_factor = Symbol("integer_factor", dimensionless, positive=True)
wavelength = Symbol("wavelength", units.length, positive=True)
string_length = Symbol("string_length", units.length, positive=True)

law = Eq(integer_factor * wavelength / 2, string_length)

# TODO: derive from boundary condition `u(L, t) = 0`


def print_law() -> str:
    return print_expression(law)


@validate_input(
    integer_factor_=integer_factor,
    string_length_=string_length,
)
@validate_output(wavelength)
def calculate_wavelength(
    integer_factor_: int,
    string_length_: Quantity,
) -> Quantity:
    result_expr = solve(law, wavelength)[0]
    result = result_expr.subs({
        integer_factor: integer_factor_,
        string_length: string_length_,
    })
    return Quantity(result)
