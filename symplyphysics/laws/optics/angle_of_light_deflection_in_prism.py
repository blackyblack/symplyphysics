from sympy import (Eq, solve)
from symplyphysics import (Quantity, Symbol, print_expression, validate_input, validate_output,
    dimensionless, angle_type)

# Description
## A prism, an optical prism, is a body made of a homogeneous material transparent to optical radiation,
## bounded by flat reflecting and refractive surfaces located at strictly defined angles to each other.
## With a small angle of incidence of the incoming beam, the angle of deflection of the beam depends only
## on the angle between the faces of the prism and the refractive index of the prism.
## There is a ray falling on the prism and a ray coming out of the prism. Let's continue these rays inside
## the prism. Then the angle of intersection of these rays, looking towards the output beam, will be called
## the angle of deviation.
## https://ru.wikipedia.org/wiki/Призма_(оптика)#:~:text=Призма%2C%20оптическая%20призма%20—%20тело%20из,определёнными%20углами%20друг%20к%20другу.

## Law is: b = a * (n - 1), where
## b - angle deviation,
## a - angle between faces,
## n - refractive index of prism.

# Conditions:
## - angle deviation of prism and angle of incidence of incoming beam are small.

angle_deviation = Symbol("angle_deviation", angle_type)

angle_faces = Symbol("angle_faces", angle_type)
refractive_index = Symbol("refractive_index", dimensionless)

law = Eq(angle_deviation, (angle_faces * (refractive_index - 1)))


def print_law() -> str:
    return print_expression(law)


@validate_input(angle_faces_=angle_faces, refractive_index_=refractive_index)
@validate_output(angle_deviation)
def calculate_angle_deviation(angle_faces_: float | Quantity, refractive_index_: float) -> Quantity:
    result_expr = solve(law, angle_deviation, dict=True)[0][angle_deviation]
    result_expr = result_expr.subs({
        angle_faces: angle_faces_,
        refractive_index: refractive_index_,
    })
    return Quantity(result_expr)
