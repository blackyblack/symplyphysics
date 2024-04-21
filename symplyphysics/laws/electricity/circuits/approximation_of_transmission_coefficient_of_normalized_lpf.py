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
## Fn(w) - approximating function of the order of n,
## w - frequency normalized to the cutoff frequency(the frequency that corresponds to the maximum bandwidth distortion).

frequency = Symbol("frequency", units.frequency)
filter_function = Function("filter_function", dimensionless)
bandwidth_distortion = Symbol("bandwidth_distortion", dimensionless)
transmission_coefficient = Symbol("transmission_coefficient", dimensionless)

law = Eq(transmission_coefficient, 1 / (1 + bandwidth_distortion**2 * filter_function(frequency)**2))


def print_law() -> str:
    return print_expression(law)


def apply_filter_function(filter_function_: Expr) -> Expr:
    applied_law = law.subs(filter_function(frequency), filter_function_)
    return applied_law


@validate_input(frequency_=frequency,
                bandwidth_distortion_=bandwidth_distortion)
@validate_output(transmission_coefficient)
def calculate_coefficient(filter_function_: Expr, frequency_: Quantity, bandwidth_distortion_: float) -> float:
    filter_function_quantity = Quantity(filter_function_)
    assert SI.get_dimension_system().equivalent_dims(filter_function_quantity.dimension,
        filter_function.dimension)
    applied_law = apply_filter_function(filter_function_)
    result_expr = applied_law.subs({
        frequency: frequency_,
        bandwidth_distortion: bandwidth_distortion_,
    })
    result = solve(result_expr, transmission_coefficient, dict=True)[0][transmission_coefficient]
    return convert_to_float(result)
