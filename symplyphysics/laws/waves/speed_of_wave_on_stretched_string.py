from sympy import Eq, sqrt
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Symbol,
    Quantity,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The speed of a wave on a stretched ideal string is set by properties of the string
## and not by properties of the wave such as frequency and amplitude.

# Law: v = sqrt(tau / mu)
## v - phase speed of wave
## tau - string tension
## mu - linear density of string

wave_speed = Symbol("phase_speed", units.velocity)
tension_force = clone_symbol(symbols.dynamics.force, "tension_force")
string_linear_density = Symbol("string_linear_density", units.mass / units.length)

law = Eq(wave_speed, sqrt(tension_force / string_linear_density))

# TODO: derive from Newton's second law


def print_law() -> str:
    return print_expression(law)


@validate_input(
    string_tension_=tension_force,
    string_linear_density_=string_linear_density,
)
@validate_output(wave_speed)
def calculate_wave_speed(
    string_tension_: Quantity,
    string_linear_density_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        tension_force: string_tension_,
        string_linear_density: string_linear_density_,
    })
    return Quantity(result)
