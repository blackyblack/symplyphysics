from sympy import Eq, solve
from symplyphysics import Symbol, units, print_expression, Quantity, \
    validate_input, validate_output

# Description
## The optical power of a lens is a value characterizing the refractive power
## of axisymmetric lenses and centered optical systems made of such lenses.

# Law D = 1 / F
# Where:
## F - focus distance of lens
## D - optical power of lens

focus_distance = Symbol("focus_distance", units.length)
optical_power = Symbol("optical_power", 1 / units.length)

law = Eq(optical_power, 1 / focus_distance)


def print_law() -> str:
    return print_expression(law)


@validate_input(focus_distance_=focus_distance)
@validate_output(optical_power)
def calculate_optical_power(focus_distance_: Quantity) -> Quantity:
    solved = solve(law, optical_power, dict=True)[0][optical_power]
    result_expr = solved.subs({
        focus_distance: focus_distance_
    })
    result = Quantity(result_expr)
    return result
