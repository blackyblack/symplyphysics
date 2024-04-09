from sympy import Eq, solve, pi, ln
from sympy.physics.units import magnetic_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)

## Description
## A coaxial waveguide is an electrical cable consisting of a central conductor and a shield arranged coaxially and separated
## by an insulating material or an air gap. It is used to transmit radio frequency electrical signals.
## The specific inductance of a coaxial waveguide depends on the radius of the outer conductor and the radius of the inner conductor,
## as well as on the relative permeability of the insulator material.

## Law is: L = (mu0 * mur / (2 * pi)) * ln(b / a), where
## L - specific inductance of coaxial waveguide,
## mu0 - magnetic constant,
## mur - relative permeability of the insulator material,
## b - radius of the outer conductor,
## a - radius of the inner conductor.

specific_inductance = Symbol("specific_inductance", units.inductance / units.length)

relative_permeability = Symbol("relative_permeability", dimensionless)
outer_radius = Symbol("outer_radius", units.length)
inner_radius = Symbol("inner_radius", units.length)

law = Eq(specific_inductance,
    (magnetic_constant * relative_permeability / (2 * pi)) * ln(outer_radius / inner_radius))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permeability_=relative_permeability,
    outer_radius_=outer_radius,
    inner_radius_=inner_radius)
@validate_output(specific_inductance)
def calculate_specific_inductance(relative_permeability_: float, outer_radius_: Quantity,
    inner_radius_: Quantity) -> Quantity:
    if outer_radius_.scale_factor <= inner_radius_.scale_factor:
        raise ValueError("The outer radius must be greater than the inner radius")
    result_velocity_expr = solve(law, specific_inductance, dict=True)[0][specific_inductance]
    result_expr = result_velocity_expr.subs({
        relative_permeability: relative_permeability_,
        outer_radius: outer_radius_,
        inner_radius: inner_radius_
    })
    return Quantity(result_expr)
