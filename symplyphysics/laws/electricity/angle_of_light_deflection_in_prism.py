from sympy import (Eq, solve)
from symplyphysics import (
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
    angle_type
)

# Description
## A prism, an optical prism, is a body made of a homogeneous material transparent to optical radiation,
## bounded by flat reflecting and refractive surfaces located at strictly defined angles to each other.
## With a small angle of incidence of the incoming beam, the angle of deflection of the beam depends only
## on the angle between the faces of the prism and the refractive index of the prism.

## Law is: b = a * (n - 1), where
## b - angle deviation,
## a - angle between faces,
## n - refractive index of prism.

angle_deviation = Symbol("angle_deviation", angle_type)

angle_faces = Symbol("angle_faces", angle_type)
refractive_index = Symbol("refractive_index", dimensionless)

law = Eq(angle_deviation, (angle_faces * (refractive_index - 1)))


def print_law() -> str:
    return print_expression(law)


@validate_input(angle_faces_=angle_faces, refractive_index_=refractive_index)
@validate_output(angle_deviation)
def calculate_angle_deviation(angle_faces_: float | Quantity, refractive_index_: float) -> float | Quantity:
    result_expr = solve(law, angle_deviation, dict=True)[0][angle_deviation]
    result_expr = result_expr.subs({
        angle_faces: angle_faces_,
        refractive_index: refractive_index_,
    })
    return Quantity(result_expr)
