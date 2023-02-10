from symplyphysics import (
    symbols, Eq, pretty, solve, units, Quantity,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Any optic lens creates image of an object. Dbstances lens-object and lens-image depend on lensoptical strength.
## This law is also called thin lens formula.

## Law: 1 / F = 1 / d + 1 / f
## Where:
## F is focus of lens,
## d is distance from lens to object,
## f is distance from lens to image.

# Conditions
## Lens is thin - it's thickness is much less than F, f and d.

focus = symbols('focus')
distance_to_object = symbols('distance_to_object')
distance_to_image = symbols('distance_to_image')

law = Eq((1 / focus), (1 / distance_to_object) + (1 / distance_to_image))

def print():
    return pretty(law, use_unicode=False)

@validate_input(object_distance_=units.length, image_distance_=units.length)
@validate_output(units.length)
def calculate_focus(object_distance_: Quantity, image_distance_: Quantity) -> Quantity:        
    result_expr = solve(law, focus, dict=True)[0][focus]
    result_expr_substituted = result_expr.subs({distance_to_object: object_distance_, distance_to_image: image_distance_})
    return expr_to_quantity(result_expr_substituted, 'focus')
