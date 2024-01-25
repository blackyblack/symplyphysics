from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output, dimensionless)

# Description
## Optical path length is the length that light needs to travel through a vacuum to create the same phase difference as it would have when traveling through a given medium. 
## If two light beams have common starting and ending points, then the difference in the optical path lengths of such rays is called the optical path difference.
## The the optical path difference corresponds to the phase shift undergone by the light emitted from two previously coherent sources when passed through mediums of different refractive indices.
## If the difference in the course of the waves is equal to an integer number of waves (i.e. an even number of half-waves), then an interference maximum is formed at the point of superposition of these waves. 
## If the difference in the course of the waves is equal to an odd number of half-waves, then an interference minimum is formed at the point of superposition of these waves.

## Law: L = l * n
## Where:
## L is optical path length
## l is geometric path length
## n is refractive index

## Conditions
## The environment should be homogeneous.


optical_path = Symbol("optical_path", units.length)
geometric_path = Symbol("geometric_path", units.length)
refractive_index = Symbol("refractive_index", dimensionless)

law = Eq(optical_path, geometric_path * refractive_index)


def print_law() -> str:
    return print_expression(law)


@validate_input(geometric_path_=geometric_path, refractive_index_=refractive_index)
@validate_output(optical_path)
def calculate_optical_path(geometric_path_: Quantity, refractive_index_: float) -> Quantity:
    result_expr = solve(law, optical_path, dict=True)[0][optical_path]
    result_optical_path = result_expr.subs({
        geometric_path: geometric_path_,
        refractive_index: refractive_index_,
    })
    return Quantity(result_optical_path)
