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
## Linear (sometimes called lateral or transverse) magnification refers to the ratio of image length to object length measured in planes that are perpendicular to the optical axis.

## Law: M = h_image / h_object
## Where:
## M is linear magnification produced by lense
## h_image is height of the image
## h_object is height of the object

## Conditions
## If the image is straight, then the value h_image is considered with a plus sign, and if it is inverted, then with a minus sign.

image_height = Symbol("image_height", units.length)
object_height = Symbol("object_height", units.length)
magnification = Symbol("magnification", dimensionless)

law = Eq(magnification, image_height / object_height)


def print_law() -> str:
    return print_expression(law)


@validate_input(image_height_=image_height, object_height_=object_height)
@validate_output(magnification)
def calculate_magnification(image_height_: Quantity, object_height_: Quantity) -> float:
    result_expr = solve(law, magnification, dict=True)[0][magnification]
    result_magnification = result_expr.subs({
        image_height: image_height_,
        object_height: object_height_,
    })
    result = Quantity(result_magnification)
    return convert_to_float(Quantity(result))
