from sympy import Eq, sin, solve
from symplyphysics import (
    angle_type,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## If the speed of a source relative to the medium exceeds the speed of sound in the medium,
## the Doppler equation no longer applies and this results in shock waves. The wavefronts
## of the waves originating from the source form a cone, namely a Mach cone. The half-angle
## of the cone is called the Mach cone angle, which is related to the Mach number of the source.
## See [this for the illustration](https://www.grc.nasa.gov/www/k-12/airplane/machang.html) of the phenomenon.

# Law: sin(theta) = 1 / M
## theta - Mach cone angle, i.e. the angle between the Mach wave wavefront (the Mach cone) and
##         the vector pointing opposite to the vector of motion
## M - Mach number of moving source

# Condition
## - M >= 1, i.e. the source moves with a speed greater than the speed of sound in the medium

mach_cone_angle = Symbol("mach_cone_angle", angle_type, positive=True)
mach_number = Symbol("mach_number", dimensionless, positive=True)

law = Eq(sin(mach_cone_angle), 1 / mach_number)


def print_law() -> str:
    return print_expression(law)


@validate_input(mach_number_=mach_number)
@validate_output(mach_cone_angle)
def calculate_mach_cone_angle(mach_number_: float) -> Quantity:
    if mach_number_ < 1:
        raise ValueError("The Mach number must be greater or equal to 1")

    result_expr = solve(law, mach_cone_angle)[1]
    result = result_expr.subs(mach_number, mach_number_)
    return Quantity(result)
