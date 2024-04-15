from sympy import Eq, solve, pi, acos, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float
)

## Description
## The relative operating bandwidth of a quarter-wave transformer depends on the reflection coefficient,
## the characteristic resistance of the transmission line and the load resistance.
## The relative operating bandwidth is the ratio of the bandwidth to the center frequency.


relative_bandwidth = Symbol("relative_bandwidth", dimensionless)

load_resistance = Symbol("load_resistance", units.impedance)
characteristic_resistance = Symbol("characteristic_resistance", units.impedance)
reflection_coefficient = Symbol("reflection_coefficient", dimensionless)

law = Eq(relative_bandwidth, 2 - (4 / pi) * acos(reflection_coefficient * 2 * sqrt(load_resistance * characteristic_resistance) / (sqrt(1 - reflection_coefficient**2) * abs(load_resistance - characteristic_resistance))))


def print_law() -> str:
    return print_expression(law)


@validate_input(load_resistance_=load_resistance,
    characteristic_resistance_=characteristic_resistance,
    reflection_coefficient_=reflection_coefficient)
@validate_output(relative_bandwidth)
def calculate_relative_bandwidth(load_resistance_: Quantity, characteristic_resistance_: Quantity,
    reflection_coefficient_: Quantity) -> float:
    result_expr = solve(law, relative_bandwidth, dict=True)[0][relative_bandwidth]
    result_expr = result_expr.subs({
        load_resistance: load_resistance_,
        characteristic_resistance: characteristic_resistance_,
        reflection_coefficient: reflection_coefficient_
    })
    return convert_to_float(Quantity(result_expr))
