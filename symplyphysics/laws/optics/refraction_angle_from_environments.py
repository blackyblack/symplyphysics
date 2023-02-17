from symplyphysics import (
    symbols, Eq, pretty, solve, units, Quantity, sin,
    validate_input, validate_output, expr_to_quantity
)

from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

# Description
## If ray of light comes from one media to another, it refracts.
## Incoming ray, refracted ray and perpendicular to media boundary are within the same flat.

## Law: n1 * sin(alpha) = n2 * sin(beta)
## Where:
## n1 is refractive index of first media,
## n2 is refractive index of second media,
## alfa is incoming angle,
## beta is refraction angle.

# Conditions
## For the sake of sanity alfa and beta are from 0 to pi/2.

refractive_index_1 = symbols('refractive_index_1')
refractive_index_2 = symbols('refractive_index_2')
incoming_angle = symbols('incoming_angle')
refraction_angle = symbols('refraction_angle')

law = Eq(refractive_index_1 * sin(incoming_angle), refractive_index_2 * sin(refraction_angle))

def print():
    return pretty(law, use_unicode=False)

@validate_input(index_1_=units.amount, incoming_angle_=angle_type, index_2_=units.amount)
@validate_output(angle_type)
def calculate_angle(index_1_: Quantity, incoming_angle_: Quantity, index_2_: Quantity) -> Quantity:        
    result_expr = solve(law, refraction_angle, dict=True)[0][refraction_angle]
    #HACK: sympy angles are always in radians
    incoming_angle_radians = incoming_angle_.scale_factor    
    angle_applied = result_expr.subs({refractive_index_1: index_1_, incoming_angle: incoming_angle_radians, refractive_index_2: index_2_})
    return expr_to_quantity(angle_applied, 'refraction_angle')
