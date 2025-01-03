from sympy import Eq, solve

from symplyphysics import Symbol, units, Quantity, validate_input, validate_output

# Description
## The optical difference in the course of two rays is the difference
## in the optical distances traversed by each of the rays

# Law: delta = L2 - L1
# Where:
# L1 - optical distance for 1st wave
# L2 - optical distance for 2nd wave
# TODO: elaborate on optical distance

# Links: Wikipedia <https://en.wikipedia.org/wiki/Optical_path_length#Optical_path_difference>

optical_distance1 = Symbol("optical_distance1", units.length)
optical_distance2 = Symbol("optical_distance2", units.length)
optical_difference_distance = Symbol("optical_difference_distance", units.length)

law = Eq(optical_difference_distance, optical_distance2 - optical_distance1)


@validate_input(optical_distance1_=optical_distance1, optical_distance2_=optical_distance2)
@validate_output(optical_difference_distance)
def calculate_optical_difference_distance(optical_distance1_: Quantity,
    optical_distance2_: Quantity) -> Quantity:
    solved = solve(law, optical_difference_distance, dict=True)[0][optical_difference_distance]
    result_expr = solved.subs({
        optical_distance1: optical_distance1_,
        optical_distance2: optical_distance2_
    })
    return Quantity(result_expr)
