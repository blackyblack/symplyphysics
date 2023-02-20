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
'''
def print():
    return pretty(law, use_unicode=False)
'''
incedence_angle_ = units.Quantity('incedence_angle')
SI.set_quantity_dimension(incedence_angle_, angle_type)
SI.set_quantity_scale_factor(incedence_angle_, 30 * units.degree)

result_expr = solve(law, refraction_angle, dict=True)[0][refraction_angle]
print(result_expr)
incoming_angle_radians = incedence_angle_.scale_factor    
print(incoming_angle_radians)
angle_applied = result_expr.subs({incedence_refractive_index: 1.003, incedence_angle: incoming_angle_radians, resulting_refractive_index: 1.333})
print(angle_applied)
# 

@validate_input(incoming_angle_=angle_type)
@validate_output(angle_type)
def calculate_angle(index_1_, incoming_angle_: Quantity, index_2_) -> Quantity:        
    result_expr = solve(law, refraction_angle, dict=True)[0][refraction_angle]
    #HACK: sympy angles are always in radians
    incoming_angle_radians = incoming_angle_.scale_factor    
    angle_applied = result_expr.subs({incedence_refractive_index: index_1_, incedence_angle: incoming_angle_radians, resulting_refractive_index: index_2_})
    return expr_to_quantity(angle_applied, 'refraction_angle')
