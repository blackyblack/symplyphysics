from symplyphysics import (
    symbols, Eq, pretty, solve, units, Quantity, sin, pi,
    validate_input, validate_output, expr_to_quantity, SI
)

from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

# Description
## If ray of light comes from one media to another, it refracts.
## Incoming ray, refracted ray and perpendicular to media boundary are within the same plane.

## Law: n1 * sin(alpha) = n2 * sin(beta)
## Where:
## n1 is refractive index of first media,
## n2 is refractive index of second media,
## alfa is incoming angle,
## beta is refraction angle.

# Conditions
## For the sake of sanity alfa and beta are from 0 to pi/2.

incedence_refractive_index = symbols('incedence_refractive_index')
resulting_refractive_index = symbols('resulting_refractive_index')
incedence_angle = symbols('incedence_angle')
refraction_angle = symbols('refraction_angle')

law = Eq(incedence_refractive_index * sin(incedence_angle), resulting_refractive_index * sin(refraction_angle))

def print():
    return pretty(law, use_unicode=False)

@validate_input(incedence_angle_=angle_type)
@validate_output(angle_type)
def calculate_refraction_angle(incedence_angle_: Quantity, incedence_refractive_index_, resulting_refractive_index_) -> Quantity:
    result_expr = solve(law, refraction_angle, dict=True)[0][refraction_angle]
    #HACK: sympy angles are always in radians
    incedence_angle_radians = incedence_angle_.scale_factor
    angle_applied = result_expr.subs({
        incedence_angle: incedence_angle_radians,
        incedence_refractive_index: incedence_refractive_index_,
        resulting_refractive_index: resulting_refractive_index_})
    #HACK: there are 2 solutions for refraction_angle: asin() and pi - asin(). We choose latter and make a sanity
    #      check for resulting angle.
    if(angle_applied > pi/2):
        angle_applied = pi - angle_applied
    resulting_angle = units.Quantity('resulting_angle')
    SI.set_quantity_dimension(resulting_angle, angle_type)
    SI.set_quantity_scale_factor(resulting_angle, angle_applied * units.radian)
    return resulting_angle
