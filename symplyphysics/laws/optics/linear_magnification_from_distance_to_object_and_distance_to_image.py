from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

# Description
## Magnification, in optics, the size of an image relative to the size of the object creating it.
## Depending on the position of the object in relation to the lens, the linear dimensions of the image change.
## If magnfication is positive, it means image formed is virtual and erect.
## If magnfication is negative, it means image formed is real and inverted.

## Law: M = f / d
## Where:
## M is linear magnification produced by lense
## f is distance from lens to image
## d is distance from lens to object

## Conditions
## If virtual image is formed, f is negative.
## If real image is formed, f is positive.
## d is always negative as object is on left side of the lens.

distance_to_object = Symbol("distance_to_object", units.length)
distance_to_image = Symbol("distance_to_image", units.length)
magnification = Symbol("magnification", dimensionless)

law = Eq(magnification, distance_to_image / distance_to_object)


def print_law() -> str:
    return print_expression(law)


@validate_input(distance_to_image_=distance_to_image, distance_to_object_=distance_to_object)
@validate_output(magnification)
def calculate_magnification(distance_to_image_: Quantity, distance_to_object_: Quantity) -> float:
    if distance_to_object_.scale_factor > 0:
        raise ValueError("The distance to the object must be non-positive.")
    result_expr = solve(law, magnification, dict=True)[0][magnification]
    result_magnification = result_expr.subs({
        distance_to_image: distance_to_image_,
        distance_to_object: distance_to_object_,
    })
    result = Quantity(result_magnification)
    return convert_to_float(result)
