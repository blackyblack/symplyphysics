from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity, sin, cos
)

# Description
## Most of cases might be represented in 2-dimensional space with two orthogonal axis - vertical Y and horizontal X. Any vector in this space (velocity, force etc) can be easily
## transformed to sum of two orthogonal vectors - projections of this vector to axis with help of angle between the vector and the horizontal axis X.

velocity  = symbols('velocity')
angle = symbols('angle')
vertical_velocity = symbols('vertical_velocity')
horizontal_velocity = symbols('horizontal_velocity')

vertical_projection = Eq(vertical_velocity, velocity * sin(angle))
horizontal_projection = Eq(horizontal_velocity, velocity, cos(angle))

def print():
    return pretty(vertical_projection, horizontal_projection, use_unicode=False)


# is there any way to make this law independent from units? we just need to check if input has the same units as output
@validate_input(object_mass_=units.mass)
@validate_output(units.force)
def calculate_horizontal_projection(generator_mass_: Quantity, object_mass_: Quantity, distance_between_mass_centers_: Quantity) -> Quantity:
    result_projection_expr = solve(law, gravity_force, dict=True)[0][gravity_force]
    result_expr = result_force_expr.subs({generator_mass: generator_mass_, object_mass: object_mass_, distance_between_mass_centers: distance_between_mass_centers_})
    return expr_to_quantity(result_expr, 'gravity_force')
