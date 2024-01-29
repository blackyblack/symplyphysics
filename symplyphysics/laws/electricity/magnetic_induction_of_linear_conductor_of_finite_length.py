from sympy import (Eq, solve, pi, cos)
from sympy.physics.units import magnetic_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, angle_type)

# Description
## Let there be a rectilinear conductor of finite length. Then its magnetic induction will depend on
## the magnitude of the current and the material. It also depends on the perpendicular distance to the
## conductor and on the angles between the lines drawn from the ends of the conductor to the point and
## the conductor.

## Law is: B = mu * mu0 * I * (cos(a1) + cos(a2)) / (4 * pi * r), where
## B - induction,
## mu - relative permeability of medium,
## mu0 - magnetic constant,
## I - current,
## a1 - first angle,
## a2 - second angle,
## r - distance.

# Conditions:
## - Conductor should be rectilinear;
## - Length of the conductor is finite.

induction = Symbol("induction", units.magnetic_density)

relative_permeability = Symbol("relative_permeability", dimensionless)
current = Symbol("current", units.current)
first_angle = Symbol("first_angle", angle_type)
second_angle = Symbol("second_angle", angle_type)
distance = Symbol("distance", units.length)

law = Eq(
    induction,
    relative_permeability * magnetic_constant * current * (cos(first_angle) + cos(second_angle)) /
    (4 * pi * distance))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permeability_=relative_permeability,
    current_=current,
    first_angle_=first_angle,
    second_angle_=second_angle,
    distance_=distance)
@validate_output(induction)
def calculate_induction(relative_permeability_: float, current_: Quantity,
    first_angle_: float | Quantity, second_angle_: float | Quantity,
    distance_: Quantity) -> Quantity:
    result_expr = solve(law, induction, dict=True)[0][induction]
    result_expr = result_expr.subs({
        relative_permeability: relative_permeability_,
        current: current_,
        first_angle: first_angle_,
        second_angle: second_angle_,
        distance: distance_
    })
    return Quantity(result_expr)
