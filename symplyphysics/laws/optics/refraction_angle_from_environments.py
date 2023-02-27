from symplyphysics import (
    symbols, Eq, pretty, solve, units, Quantity, sin, pi,
    validate_input, validate_output, SI
)

from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

# Description
## If ray of light comes from one media to another, it refracts.
## Incoming ray, refracted ray and perpendicular to media boundary are within the same plane.

## Law: n1 * sin(alpha) = n2 * sin(beta)
## Where:
## n1 is refractive index of first media.
## n2 is refractive index of second media.
## alfa is incoming angle. Angle is measured between normal vector to a media boundary plane directed outside of refracting media and incoming ray of light.
## beta is refraction angle. Angle is measured between normal vector to a media boundary plane directed inside refracting media and refracted ray of light.

# Conditions
## - alfa is in [-pi/2, pi/2] range, meaning that incoming ray comes from the outer media.
## - beta is in [-pi/2, pi/2] range, meaning that refracted ray stays in the refracting media.
## - light is monochromic, as refactive index depends on the light frequency.
## - refracting media is uniform, so refracting index does not change over ray path.

incedence_refractive_index = symbols('incedence_refractive_index')
resulting_refractive_index = symbols('resulting_refractive_index')
incedence_angle = symbols('incedence_angle')
refraction_angle = symbols('refraction_angle')

law = Eq(incedence_refractive_index * sin(incedence_angle), resulting_refractive_index * sin(refraction_angle))

def print():
    return pretty(law, use_unicode=False)

@validate_input(incedence_angle_=angle_type)
@validate_output(angle_type)
def calculate_refraction_angle(incedence_angle_: Quantity, incedence_refractive_index_: float, resulting_refractive_index_: float) -> Quantity:
    #HACK: sympy angles are always in radians
    incedence_angle_radians = incedence_angle_.scale_factor
    # Check for boundary conditions
    assert incedence_angle_radians <= pi/2
    assert incedence_angle_radians >= -pi/2
    solutions = solve(law, refraction_angle, dict=True)
    result_expr = solutions[0][refraction_angle]
    angle_applied = result_expr.subs({
        incedence_angle: incedence_angle_radians,
        incedence_refractive_index: incedence_refractive_index_,
        resulting_refractive_index: resulting_refractive_index_})
    #HACK: there are 2 solutions for refraction_angle: pi - asin() and asin(). We choose former and switch to latter
    #      if resulting angle does not pass boundary conditions.
    if(angle_applied > pi/2 or angle_applied < -pi/2):
        result_expr = solutions[1][refraction_angle]
        angle_applied = result_expr.subs({
            incedence_angle: incedence_angle_radians,
            incedence_refractive_index: incedence_refractive_index_,
            resulting_refractive_index: resulting_refractive_index_})

    # Check for boundary conditions
    assert angle_applied <= pi/2
    assert angle_applied >= -pi/2
    resulting_angle = units.Quantity('resulting_angle')
    SI.set_quantity_dimension(resulting_angle, angle_type)
    SI.set_quantity_scale_factor(resulting_angle, angle_applied * units.radian)
    return resulting_angle
