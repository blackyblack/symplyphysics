from sympy import Eq, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)
from symplyphysics.core.convert import convert_to_dimensionless

# Description
# Froude number characterizes the ratio between the force of inertia and
# the external force, in the field of which the motion occurs, acting on
# an elementary volume of liquid or gas. There is a characteristic length in
# the formula. The characteristic length is the dimension
# that defines the length scale of a physical system. A characteristic length
# is usually the volume of a system divided by its surface: L = V / A,
# where V is the volume of the body, and A is the cross-sectional area.
# For example, it is used to calculate flow through circular and non-circular
# tubes in order to examine flow conditions. D = 4 * A / p, where
# D is characteristic diameter, A is the cross-sectional are, p is wetted perimeter.
# Law: Fr = u / sqrt(g * L), where
# u is velocity,
# g is acceleration due to gravity,
# L is a characteristic length,
# Fr is Froude number.

velocity = Symbol("velocity", units.velocity)
characteristic_length = Symbol("characteristic_length", units.length)
froude_number = Symbol("froude_number", dimensionless)

law = Eq(froude_number, velocity / sqrt(units.acceleration_due_to_gravity * characteristic_length))


def print_law() -> str:
    return print_expression(law)


@validate_input(velocity_=velocity, characteristic_length_=characteristic_length)
@validate_output(froude_number)
def calculate_froude_number(velocity_: Quantity, characteristic_length_: Quantity) -> float:
    result_expr = solve(law, froude_number, dict=True)[0][froude_number]
    result_applied = result_expr.subs({
        velocity: velocity_,
        characteristic_length: characteristic_length_
    })
    result = Quantity(result_applied)
    return convert_to_dimensionless(Quantity(result))
