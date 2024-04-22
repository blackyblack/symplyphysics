from math import ceil
from sympy import (
    Expr,
    Eq,
    solve,
)
from symplyphysics import (
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
## The band-stop distortion distortion determines the level of distortion in the band-stop.

# Law: F(fs, n) = es / ep
## Where:
## F(f, n) - approximating function of the order of n,
## es - band-stop distortion,
## ep - bandwidth distortion,
## fs - frequency of the band-stop.

band_stop_frequency = Symbol("band_stop_frequency", units.frequency)
filter_order = Symbol("filter_order", dimensionless)
filter_function = Function("filter_function", dimensionless)
bandwidth_distortion = Symbol("bandwidth_distortion", dimensionless)
band_stop_distortion = Symbol("band_stop_distortion", dimensionless)

law = Eq(filter_function(band_stop_frequency, filter_order), band_stop_distortion / bandwidth_distortion)


def print_law() -> str:
    return print_expression(law)


def apply_filter_function(filter_function_: Expr) -> Expr:
    applied_law = law.subs(filter_function(band_stop_frequency, filter_order), filter_function_)
    return applied_law


@validate_input(band_stop_frequency_=band_stop_frequency,
                bandwidth_distortion_=bandwidth_distortion,
                band_stop_distortion_=band_stop_distortion)
@validate_output(filter_order)
def calculate_order(filter_function_: Expr, band_stop_frequency_: Quantity, bandwidth_distortion_: float,
                          band_stop_distortion_: float) -> int:
    applied_law = apply_filter_function(filter_function_)
    result_expr = applied_law.subs({
        band_stop_frequency: band_stop_frequency_,
        bandwidth_distortion: bandwidth_distortion_,
        band_stop_distortion: band_stop_distortion_,
    })
    result = solve(result_expr, filter_order, dict=True)[0][filter_order]
    return ceil(convert_to_float(result))
