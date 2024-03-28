from sympy import (Eq, solve, S)
from sympy.physics.units import magnetic_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, convert_to, dimensionless, angle_type)


# Description
## A rectangular resonator consists of metal walls and a material filling it.
## In the case when the resonator is empty, its quality factor depends only on the losses in the metal walls of the resonator.
## The quality factor shows the ratio of the energy stored in the resonator to the energy loss during one oscillation period.

## Law is: Q = w0 * mu0 * mur * b * a * l * (a^2 + l^2) / (2 * Rs * (a^3 * (l + 2 * b) + l^3 * (a + 2 * b))), where
## Q - quality factor of the resonator,
## mur - relative permeability of the resonator walls,
## Rs1 - surface resistance of the resonator walls,
## a - height of the resonator,
## b - width of the resonator,
## l - length of the resonator.

# Conditions:
# - for transverse electric waves with indices m, n, p = 1, 0, 1.

quality_factor = Symbol("quality_factor", dimensionless)

resonant_frequency = Symbol("resonant_frequency", angle_type / units.time)
relative_permeability = Symbol("relative_permeability", dimensionless)
surface_resistance = Symbol("surface_resistance", units.impedance)
resonator_width = Symbol("resonator_width", units.length)
resonator_height = Symbol("resonator_height", units.length)
resonator_length = Symbol("resonator_length", units.length)

law = Eq(quality_factor, resonant_frequency * magnetic_constant * relative_permeability * resonator_width * resonator_height * resonator_length * (resonator_height**2 + resonator_length**2)
         / (2 * surface_resistance * (resonator_height**3 * (resonator_length + 2 * resonator_width) + resonator_length**3 * (resonator_height + 2 * resonator_width))))


def print_law() -> str:
    return print_expression(law)


@validate_input(resonant_frequency_=resonant_frequency,
    relative_permeability_=relative_permeability,
    surface_resistance_=surface_resistance,
    resonator_dimensions_=units.length)
@validate_output(quality_factor)
def calculate_quality_factor(resonant_frequency_: Quantity, relative_permeability_: float, surface_resistance_: Quantity,
    resonator_dimensions_: tuple[Quantity, Quantity, Quantity]) -> float:
    resonator_width_, resonator_height_, resonator_length_ = resonator_dimensions_
    result_expr = solve(law, quality_factor, dict=True)[0][quality_factor]
    result_expr = result_expr.subs({
        resonant_frequency: resonant_frequency_,
        relative_permeability: relative_permeability_,
        surface_resistance: surface_resistance_,
        resonator_width: resonator_width_,
        resonator_height: resonator_height_,
        resonator_length: resonator_length_
    })
    return float(convert_to(Quantity(result_expr), S.One).evalf())
