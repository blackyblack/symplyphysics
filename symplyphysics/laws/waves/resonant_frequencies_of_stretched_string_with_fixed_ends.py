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
## For a string with fixed ends there is only a limited set of frequencies at which standing waves will
## occur on it. Each possible frequency is a resonant frequency, and the corresponding wave pattern is
## an oscillation mode. The oscillation mode corresponding to n = 1 is called the fundamental mode or
## the first harmonic; the mode corresponding to n = 2 is the second harmonic; and so on.

# Law: f = n * v / (2 * L)
## f - resonant frequency of n-th harmonic
## n - integer (harmonic number)
## v - phase velocity of wave
## L - length of string

resonant_frequency = Symbol("resonant_frequency", 1 / units.time, positive=True)
harmonic_number = Symbol("harmonic_number", dimensionless, integer=True, positive=True)
phase_velocity = Symbol("phase_velocity", units.velocity, positive=True)
string_length = Symbol("string_length", units.length, positive=True)

law = Eq(resonant_frequency, harmonic_number * phase_velocity / (2 * string_length))

# TODO: Derive from the condition ofa  wave along a stretched string with fixed ends


def print_law() -> str:
    return print_expression(law)


@validate_input(
    harmonic_number_=harmonic_number,
    phase_velocity_=phase_velocity,
    string_length_=string_length,
)
@validate_output(resonant_frequency)
def calculate_resonant_frequency(
    harmonic_number_: int,
    phase_velocity_: Quantity,
    string_length_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        harmonic_number: harmonic_number_,
        phase_velocity: phase_velocity_,
        string_length: string_length_,
    })
    return Quantity(result)
