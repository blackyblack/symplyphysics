from sympy import (Eq, solve, S)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output, dimensionless, convert_to)

# Description
## Magnification, in optics, the size of an image relative to the size of the object creating it.
## Depending on the position of the object in relation to the lens, the linear dimensions of the image change.

## Law: M = f / d
## Where:
## M is linear magnification produced by lense
## f is distance from lens to image.
## d is distance from lens to object,

## Conditions
## f > 0 – for real images; f < 0 – for imaginary images. The "–" sign is placed before d only when a converging beam of rays falls on the lens.
## Then they are mentally extended to the intersection behind the lens, an imaginary light source is placed there, and the distance d is determined for it.


distance_to_object = Symbol("distance_to_object", units.length)
distance_to_image = Symbol("distance_to_image", units.length)
magnification = Symbol("magnification", dimensionless)

law = Eq(magnification, distance_to_image / distance_to_object)


def print_law() -> str:
    return print_expression(law)


@validate_input(distance_to_image_=distance_to_image, distance_to_object_=distance_to_object)
@validate_output(magnification)
def calculate_magnification(distance_to_image_ : Quantity, distance_to_object_: Quantity) -> float:
    result_expr = solve(law, magnification, dict=True)[0][magnification]
    result_magnification = result_expr.subs({
        distance_to_image: distance_to_image_,
        distance_to_object: distance_to_object_,
    })
    result = Quantity(result_magnification)
    return float(convert_to(result, S.One).evalf())
