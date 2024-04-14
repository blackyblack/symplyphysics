from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output
)

# Description
## A multistage resistance transformer consists of several sections. The length of each section depends on the wavelength at
## the upper operating frequency and the wavelength at the lower operating frequency.

## Law is: L = L1 * L2 / (2 * (L1 + L2)), where
## L - length of section,
## L1 - the wavelength at the upper frequency,
## L2 - the wavelength at the lower frequency.

length_of_section = Symbol("length_of_section", units.length)

wavelength_upper_frequency = Symbol("wavelength_upper_frequency", units.length)
wavelength_lower_frequency = Symbol("wavelength_lower_frequency", units.length)

law = Eq(length_of_section, wavelength_upper_frequency * wavelength_lower_frequency / (2 * (wavelength_upper_frequency + wavelength_lower_frequency)))


def print_law() -> str:
    return print_expression(law)


@validate_input(wavelength_upper_frequency_=wavelength_upper_frequency, wavelength_lower_frequency_=wavelength_lower_frequency)
@validate_output(length_of_section)
def calculate_length_of_section(wavelength_upper_frequency_: Quantity,
    wavelength_lower_frequency_: Quantity) -> Quantity:
    result_expr = solve(law, length_of_section, dict=True)[0][length_of_section]
    result_expr = result_expr.subs({
        wavelength_upper_frequency: wavelength_upper_frequency_,
        wavelength_lower_frequency: wavelength_lower_frequency_,
    })
    return Quantity(result_expr)
