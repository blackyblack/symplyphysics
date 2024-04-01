from sympy import Eq, solve, pi, sqrt, ln
from sympy.physics.units import electric_constant, magnetic_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)

## Description
## A coaxial waveguide is an electrical cable consisting of a central conductor and a shield arranged coaxially and separated
## by an insulating material or an air gap. It is used to transmit radio frequency electrical signals.
## The wave resistance of a coaxial waveguide depends on the radius of the outer conductor and the radius of the inner conductor,
## as well as on the relative permittivity and the relative permeability of the insulator material.

## Law is: Z = (1 / (2 * pi)) * sqrt(mu0 * mur / (e0 * er)) * ln(b / a), where
## Z - wave resistance of coaxial waveguide,
## e0 - electric constant,
## er - relative permittivity of insulating material,
## mu0 - magnetic constant,
## mur - relative permeability of the insulator material,
## b - radius of the outer conductor,
## a - radius of the inner conductor.

wave_resistance = Symbol("wave_resistance", units.impedance)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
relative_permeability = Symbol("relative_permeability", dimensionless)
outer_radius = Symbol("outer_radius", units.length)
inner_radius = Symbol("inner_radius", units.length)

law = Eq(wave_resistance, (1 / (2 * pi)) * sqrt(magnetic_constant * relative_permeability /
    (electric_constant * relative_permittivity)) * ln(outer_radius / inner_radius))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability,
    outer_radius_=outer_radius,
    inner_radius_=inner_radius)
@validate_output(wave_resistance)
def calculate_wave_resistance(relative_permittivity_: float, relative_permeability_: float,
    outer_radius_: Quantity, inner_radius_: Quantity) -> Quantity:
    if outer_radius_.scale_factor <= inner_radius_.scale_factor:
        raise ValueError("The outer radius must be greater than the inner radius")
    result_velocity_expr = solve(law, wave_resistance, dict=True)[0][wave_resistance]
    result_expr = result_velocity_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
        outer_radius: outer_radius_,
        inner_radius: inner_radius_
    })
    return Quantity(result_expr)
