from sympy import (Eq, solve, symbols)
from symplyphysics import print_expression

# Description
## Average value of the cosine of the angle in the lab system at which neutrons are scattered in the medium.

## It can be calculated for most of the neutron energies as:
## μ = 2 / (3 * A)
## Where:
## A is the mass number of target nucleus.
## μ is average value of the cosine of the angle in the lab system at which neutrons are scattered in the medium.

target_nucleus_mass_number = symbols("target_nucleus_mass_number")
average_scattering_angle_cosine = symbols("average_scattering_angle_cosine")

law = Eq(average_scattering_angle_cosine, 2 / (3 * target_nucleus_mass_number))


def print() -> str:
    return print_expression(law)


def calculate_average_scattering_angle_cosine(target_nucleus_mass_number_: int) -> float:
    result_angle_cosine_expr = solve(law, average_scattering_angle_cosine,
        dict=True)[0][average_scattering_angle_cosine]
    result_expr = result_angle_cosine_expr.subs(target_nucleus_mass_number,
        target_nucleus_mass_number_)
    return result_expr.evalf()
