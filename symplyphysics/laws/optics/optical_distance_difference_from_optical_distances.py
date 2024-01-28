from sympy import Eq, solve

from symplyphysics import Symbol, units, print_expression, Quantity, \
    validate_input, validate_output

# Description
## The optical difference in the course of two rays is the difference
## in the optical distances traversed by each of the rays

# Law: delta = L2 - L1
# Where:
# L1 - optical distance for 1st wave
# L2 - optical distance for 2nd wave

optical_distance_wave_1 = Symbol("", units.length)
optical_distance_wave_2 = Symbol("", units.length)
optical_difference_distance = Symbol("", units.length)

law = Eq(optical_difference_distance, optical_distance_wave_2 - optical_distance_wave_1)


def print_law() -> str:
    return print_expression(law)


@validate_input(optical_distance_wave_1_=optical_distance_wave_1, optical_distance_wave_2_=optical_distance_wave_2)
@validate_output(optical_difference_distance)
def calculate_optical_difference_distance(optical_distance_wave_1_: Quantity, optical_distance_wave_2_: Quantity) -> Quantity:
    solved = solve(law, optical_difference_distance, dict=True)[0][optical_difference_distance]
    result_expr = solved.subs({
        optical_distance_wave_1: optical_distance_wave_1_,
        optical_distance_wave_2: optical_distance_wave_2_
    })
    result = Quantity(result_expr)
    return result

