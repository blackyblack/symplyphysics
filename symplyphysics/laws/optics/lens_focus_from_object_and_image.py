from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)

# Description
## Any optic lens creates image of an object. Distances lens-object and lens-image depend on lens optical strength.
## This law is also called thin lens formula.

## Law: 1 / F = 1 / d + 1 / f
## Where:
## F is focus distance of lens,
## d is distance from lens to object,
## f is distance from lens to image.

# Conditions
## Lens is thin - it's thickness is much less than F, f and d.

# Links: Wikipedia <https://en.wikipedia.org/wiki/Lens#Lens_equation>

focus_distance = Symbol("focus_distance", units.length)
distance_to_object = Symbol("distance_to_object", units.length)
distance_to_image = Symbol("distance_to_image", units.length)

law = Eq((1 / focus_distance), (1 / distance_to_object) + (1 / distance_to_image))


@validate_input(object_distance_=distance_to_object, image_distance_=distance_to_image)
@validate_output(focus_distance)
def calculate_focus(object_distance_: Quantity, image_distance_: Quantity) -> Quantity:
    result_expr = solve(law, focus_distance, dict=True)[0][focus_distance]
    focus_applied = result_expr.subs({
        distance_to_object: object_distance_,
        distance_to_image: image_distance_
    })
    return Quantity(focus_applied)
