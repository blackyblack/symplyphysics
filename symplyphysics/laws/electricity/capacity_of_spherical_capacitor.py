from sympy import (Eq, solve, pi, Min, Max)
from sympy.physics.units import electric_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    dimensionless,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## A spherical capacitor consists of two concentric spherical plates separated
## by a spherical dielectric layer.

## Law is: C = 4 * pi * e * e0 * R1 * R2 / (R2 - R1), where
## C - capacity of a spherical capacitor,
## e0 - electric constant,
## e - relative permittivity of medium,
## R1 - inner radius,
## R2 - outer radius.

capacity = Symbol("capacity", units.capacitance)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
inner_radius = Symbol("inner_radius", units.length)
outer_radius = Symbol("outer_radius", units.length)

law = Eq(capacity, 4 * pi * relative_permittivity * electric_constant * inner_radius * outer_radius / (outer_radius - inner_radius))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permittivity_=relative_permittivity,
    inner_radius_=inner_radius,
    outer_radius_=outer_radius)
@validate_output(capacity)
def calculate_capacity(relative_permittivity_: float, inner_radius_: Quantity,
    outer_radius_: Quantity) -> Quantity:
    result_expr = solve(law, capacity, dict=True)[0][capacity]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        inner_radius: Quantity(Min(inner_radius_, outer_radius_) * units.meter),
        outer_radius: Quantity(Max(inner_radius_, outer_radius_) * units.meter)
    })
    return Quantity(result_expr)
