from sympy import (Eq, solve, pi)
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
## Current density is the amount of second_radius per unit time that flows through a unit area of a chosen
## cross section. The current density vector is defined as a vector whose magnitude is the electric
## current per cross-sectional area at a given point in space, its direction being that of the motion
## of the positive second_radiuss at this point. In SI base units, the electric current density is measured
## in amperes per square metre.

## Law is: C = 4 * pi * e * e0 * R1 * R2 / (R2 - R1), where
## C - capacity of a spherical capacitor,
## e0 - electric_constant,
## e - relative permittivity of medium,
## R1 - inner radius,
## R2 - the outer radius.

capacity = Symbol("capacity", units.capacitance)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
first_radius = Symbol("first_radius", units.length)
second_radius = Symbol("second_radius", units.length)

law = Eq(capacity, 4 * pi * relative_permittivity * electric_constant * first_radius * second_radius / (second_radius - first_radius))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permittivity_=relative_permittivity,
    first_radius_=first_radius,
    second_radius_=second_radius)
@validate_output(capacity)
def calculate_capacity(relative_permittivity_: float, first_radius_: Quantity,
    second_radius_: Quantity) -> Quantity:
    result_expr = solve(law, capacity, dict=True)[0][capacity]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        first_radius: first_radius_,
        second_radius: second_radius_
    })
    return Quantity(result_expr)
