from symplyphysics import (
    symbols, Eq, pretty, Quantity, units, solve,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Centripetal acceleration is type of acceleration which is perpendicular to velocity vector and directed to the curve center.
## Centripetal acceleration only changes the velocity vector direction and not the velocity vector length.
## Definition: an = V**2 / R, where
## an is centripetal acceleration (aka normal acceleration),
## V is linear velocity,
## R is curve radius.

centripetal_acceleration = symbols("centripetal_acceleration")
linear_velocity = symbols("linear_velocity")
curve_radius = symbols("curve_radius")

definition = Eq(centripetal_acceleration, linear_velocity**2 / curve_radius)

def print():
    return pretty(definition, use_unicode=False)

@validate_input(linear_velocity_=units.velocity, curve_radius_=units.length)
@validate_output(units.acceleration)
def calculate_acceleration(linear_velocity_: Quantity, curve_radius_: Quantity) -> Quantity:        
    solved = solve(definition, centripetal_acceleration, dict=True)[0][centripetal_acceleration]
    result_expr = solved.subs({
        linear_velocity:linear_velocity_,
        curve_radius:curve_radius_
    })
    return expr_to_quantity(result_expr, "acceleration")
