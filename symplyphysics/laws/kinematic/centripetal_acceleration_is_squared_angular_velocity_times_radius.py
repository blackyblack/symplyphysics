from sympy import Eq, solve
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    print_expression,
    angle_type,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic import (
    linear_velocity_from_angular_velocity_and_radius as velocities_law,
    centripetal_acceleration_is_squared_velocity_by_radius as centripetal_law,
)

# Description
## Centripetal acceleration is defined as the change in velocity tangential to the velocity vector.
## See more at [its definition through linear velocity](./centripetal_acceleration_is_squared_velocity_by_radius.py)

# Law: a_n = w**2 * r
## a_n - centripetal (or normal) acceleration
## w - angular velocity
## r - curve radius

centripetal_acceleration = clone_symbol(symbols.kinematic.acceleration, "centripetal_acceleration")
angular_velocity = Symbol("angular_velocity", angle_type / units.time)
curve_radius = Symbol("curve_radius", units.length)

law = Eq(centripetal_acceleration, angular_velocity**2 * curve_radius)

# Derive law from expression for linear velocity in circular motion

centripetal_acceleration_derived = centripetal_law.law.rhs.subs(centripetal_law.curve_radius,
    curve_radius)

velocities_law_sub = velocities_law.law.subs({
    velocities_law.linear_velocity: centripetal_law.linear_velocity,
    velocities_law.angular_velocity: angular_velocity,
    velocities_law.curve_radius: curve_radius,
})

centripetal_acceleration_derived = solve([
    Eq(centripetal_acceleration, centripetal_acceleration_derived),
    velocities_law_sub,
], (centripetal_acceleration, centripetal_law.linear_velocity),
    dict=True)[0][centripetal_acceleration]

assert expr_equals(law.rhs, centripetal_acceleration_derived)


def print_law() -> str:
    return print_expression(law)


@validate_input(angular_velocity_=angular_velocity, curve_radius_=curve_radius)
@validate_output(centripetal_acceleration)
def calculate_centripetal_acceleration(angular_velocity_: Quantity,
    curve_radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        angular_velocity: angular_velocity_,
        curve_radius: curve_radius_,
    })
    return Quantity(result)
