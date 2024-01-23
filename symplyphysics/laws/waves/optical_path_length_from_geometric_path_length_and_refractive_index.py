from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output, dimensionless)

# Description
## The optical path length L between points A and B of a transparent medium is the distance over which light (optical radiation) has spread it would be in a vacuum for the same time it takes from A to B in the medium.


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
