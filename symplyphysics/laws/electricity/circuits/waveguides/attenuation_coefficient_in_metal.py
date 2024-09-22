from sympy import (
    Eq,
    solve,
    pi,
    sqrt,
    ln,
)
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
## A coaxial waveguide is an electrical cable consisting of a central conductor and a shield arranged coaxially and separated
## by an insulating material or an air gap. It is used to transmit radio frequency electrical signals.
## The specific resistance of a coaxial waveguide depends on the diameter of the outer conductor and the diameter of the inner conductor,
## as well as on the relative permeability and the relative permittivity of the insulator material, the surface resistance of the outer conductor
## and the surface resistance of the inner conductor.
## The attenuation coefficient shows how many times the transmitted signal weakens per unit length of the coaxial waveguide.

## Law is: am = sqrt(er / mur) * ((Rs1 / d) + (Rs2 / D)) / (420 * pi * ln(D / d)), where
## am - attenuation coefficient in metal,
## er - relative permittivity of insulating material,
## mur - relative permeability of the insulator material,
## Rs1 - surface resistance of the inner conductor,
## Rs2 - surface resistance of the outer conductor,
## D - diameter of the outer conductor,
## d - diameter of the inner conductor.

attenuation_coefficient = Symbol("attenuation_coefficient", 1 / units.length)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
relative_permeability = Symbol("relative_permeability", dimensionless)
surface_resistance_outer = Symbol("surface_resistance_outer", units.impedance)
surface_resistance_inner = Symbol("surface_resistance_inner", units.impedance)
outer_diameter = Symbol("outer_diameter", units.length)
inner_diameter = Symbol("inner_diameter", units.length)

resistance = Quantity(420 * units.ohm)

law = Eq(
    attenuation_coefficient,
    sqrt(relative_permittivity / relative_permeability) *
    ((surface_resistance_inner / inner_diameter) + (surface_resistance_outer / outer_diameter)) /
    (resistance * pi * ln(outer_diameter / inner_diameter)))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability,
    surface_resistance_outer_=surface_resistance_outer,
    surface_resistance_inner_=surface_resistance_inner,
    outer_diameter_=outer_diameter,
    inner_diameter_=inner_diameter)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(relative_permittivity_: float, relative_permeability_: float,
    surface_resistance_outer_: Quantity, surface_resistance_inner_: Quantity,
    outer_diameter_: Quantity, inner_diameter_: Quantity) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    if outer_diameter_.scale_factor <= inner_diameter_.scale_factor:
        raise ValueError("The outer diameter must be greater than the inner diameter")
    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
        surface_resistance_outer: surface_resistance_outer_,
        surface_resistance_inner: surface_resistance_inner_,
        outer_diameter: outer_diameter_,
        inner_diameter: inner_diameter_
    })
    return Quantity(result_expr)
