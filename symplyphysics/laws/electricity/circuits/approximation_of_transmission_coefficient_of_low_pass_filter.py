from sympy import (
    Expr,
    Eq,
    solve,
)
from symplyphysics import (
    SI,
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float
)

# Description
## The approximation of the power transmission coefficient of a normalized low-pass filter is given by approximating functions of the order of n.
## Bandwidth distortion determines the maximum distortion in the bandwidth. In other words, to determine the level of ripples in the bandwidth.

# Law: H = 1 / (1 + e^2 * Fn(w)^2)
## Where:
## H - the transfer coefficient of the filter,
## e - bandwidth distortion,
## Fn(w) - approximating function of the order of n of the transfer coefficient,
## w - frequency.

frequency = Symbol("frequency", units.frequency)
filter_function = Function("filter_function", dimensionless)
bandwidth_distortion = Symbol("bandwidth_distortion", dimensionless)
transmission_coefficient = Symbol("transmission_coefficient", dimensionless)

law = Eq(transmission_coefficient, 1 / (1 + bandwidth_distortion**2 * filter_function(frequency)**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(filter_function_=filter_function,
                bandwidth_distortion_=bandwidth_distortion)
@validate_output(transmission_coefficient)
def calculate_coefficient(filter_function_: float, bandwidth_distortion_: float) -> float:
    result_expr = law.subs({
        filter_function(frequency): filter_function_,
        bandwidth_distortion: bandwidth_distortion_,
    })
    result = solve(result_expr, transmission_coefficient, dict=True)[0][transmission_coefficient]
    return convert_to_float(result)
