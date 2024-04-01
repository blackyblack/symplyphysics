from sympy import (
    Eq,
    solve,
    sqrt,
)
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)

# Description
## A rectangular resonator consists of metal walls and a material filling it.
## The indices show the number of half-waves that fit along the width, height and length of the resonator, respectively.

## Law is: f = c * sqrt((m / a)^2 + (n / b)^2 + (p / l)^2) / (2 * sqrt(er * mur)), where
## f - resonant frequency of resonator,
## c - speed of light,
## m - first index,
## n - second index,
## p - third index,
## a - width of the resonator,
## b - height of the resonator,
## l - length of the resonator,
## er - relative permittivity of the material filling the resonator,
## mur - relative permeability of the material filling the resonator.

resonant_frequency = Symbol("resonant_frequency", 1 / units.time)

first_index = Symbol("first_index", dimensionless)
second_index = Symbol("second_index", dimensionless)
third_index = Symbol("third_index", dimensionless)
resonator_width = Symbol("resonator_width", units.length)
resonator_height = Symbol("resonator_height", units.length)
resonator_length = Symbol("resonator_length", units.length)
relative_permittivity = Symbol("relative_permittivity", dimensionless)
relative_permeability = Symbol("relative_permeability", dimensionless)

law = Eq(
    resonant_frequency,
    speed_of_light * sqrt((first_index / resonator_width)**2 +
    (second_index / resonator_height)**2 + (third_index / resonator_length)**2) /
    (2 * sqrt(relative_permittivity * relative_permeability)))


def print_law() -> str:
    return print_expression(law)


@validate_input(indexes_=dimensionless,
    resonator_dimensions_=units.length,
    relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability)
@validate_output(resonant_frequency)
def calculate_resonant_frequency(indexes_: tuple[float, float,
    float], resonator_dimensions_: tuple[Quantity, Quantity, Quantity],
    relative_permittivity_: Quantity, relative_permeability_: Quantity) -> Quantity:
    first_index_, second_index_, third_index_ = indexes_
    resonator_width_, resonator_height_, resonator_length_ = resonator_dimensions_
    result_expr = solve(law, resonant_frequency, dict=True)[0][resonant_frequency]
    result_expr = result_expr.subs({
        first_index: first_index_,
        second_index: second_index_,
        third_index: third_index_,
        resonator_width: resonator_width_,
        resonator_height: resonator_height_,
        resonator_length: resonator_length_,
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
    })
    return Quantity(result_expr)
