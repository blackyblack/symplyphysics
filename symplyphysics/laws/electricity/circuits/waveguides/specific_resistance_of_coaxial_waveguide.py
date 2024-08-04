from sympy import (
    Eq,
    solve,
    pi,
    sqrt,
)
from sympy.physics.units import magnetic_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
    angle_type,
)

# Description
## A coaxial waveguide is an electrical cable consisting of a central conductor and a shield arranged coaxially and separated
## by an insulating material or an air gap. It is used to transmit radio frequency electrical signals.
## The specific resistance of a coaxial waveguide depends on the radius of the outer conductor and the radius of the inner conductor,
## as well as on the relative permeability of the insulator material, frequency of signal and specific conductivity of conductor.

## Law is: R = (1 / (2 * pi)) * sqrt(w * mu0 * mur / (2 * sig)) * (1 / a - 1 / b), where
## R - specific resistance of coaxial waveguide,
## w - angular frequency of signal,
## sig - specific conductivity of conductor,
## mu0 - magnetic constant,
## mur - relative permeability of the insulator material,
## b - radius of the outer conductor,
## a - radius of the inner conductor.

specific_resistance = Symbol("specific_resistance", units.impedance / units.length)

relative_permeability = Symbol("relative_permeability", dimensionless)
angular_frequency = Symbol("angular_frequency", angle_type / units.time)
specific_conductivity = Symbol("specific_conductivity", units.conductance / units.length)
outer_radius = Symbol("outer_radius", units.length)
inner_radius = Symbol("inner_radius", units.length)

law = Eq(specific_resistance,
    (1 / (2 * pi)) * sqrt(angular_frequency * magnetic_constant * relative_permeability /
    (2 * specific_conductivity)) * (1 / inner_radius - 1 / outer_radius))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permeability_=relative_permeability,
    angular_frequency_=angular_frequency,
    specific_conductivity_=specific_conductivity,
    outer_radius_=outer_radius,
    inner_radius_=inner_radius)
@validate_output(specific_resistance)
def calculate_specific_resistance(relative_permeability_: float, angular_frequency_: Quantity,
    specific_conductivity_: Quantity, outer_radius_: Quantity, inner_radius_: Quantity) -> Quantity:
    if outer_radius_.scale_factor <= inner_radius_.scale_factor:
        raise ValueError("The outer radius must be greater than the inner radius")
    result_expr = solve(law, specific_resistance, dict=True)[0][specific_resistance]
    result_expr = result_expr.subs({
        relative_permeability: relative_permeability_,
        angular_frequency: angular_frequency_,
        specific_conductivity: specific_conductivity_,
        outer_radius: outer_radius_,
        inner_radius: inner_radius_
    })
    return Quantity(result_expr)
