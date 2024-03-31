from sympy import Eq, solve, pi, ln
from sympy.physics.units import electric_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)

## Description
## A coaxial waveguide is an electrical cable consisting of a central conductor and a shield arranged coaxially and separated
## by an insulating material or an air gap. It is used to transmit radio frequency electrical signals.
## The specific capacitance of a coaxial waveguide depends on the radius of the outer conductor and the radius of the inner conductor,
## as well as on the relative permittivity of the insulator material.

## Law is: C = (2 * pi * e0 * er ) / ln(b / a), where
## C - specific capacitance of coaxial waveguide,
## e0 - electric constant,
## er - relative permittivity of insulating material,
## b - radius of the outer conductor,
## a - radius of the inner conductor.

specific_capacitance = Symbol("specific_capacitance", units.capacitance / units.length)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
outer_radius = Symbol("outer_radius", units.length)
inner_radius = Symbol("inner_radius", units.length)

law = Eq(specific_capacitance, (2 * pi * electric_constant * relative_permittivity) / ln(outer_radius / inner_radius))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permittivity_=relative_permittivity, outer_radius_=outer_radius, inner_radius_=inner_radius)
@validate_output(specific_capacitance)
def calculate_specific_capacitance(relative_permittivity_: float, outer_radius_: Quantity,
    inner_radius_: Quantity) -> Quantity:
    if outer_radius_.scale_factor <= inner_radius_.scale_factor:
        raise ValueError("The outer radius must be greater than the inner radius")
    result_velocity_expr = solve(law, specific_capacitance, dict=True)[0][specific_capacitance]
    result_expr = result_velocity_expr.subs({
        relative_permittivity: relative_permittivity_,
        outer_radius: outer_radius_,
        inner_radius: inner_radius_
    })
    return Quantity(result_expr)
