from sympy import Eq, Rational, solve, S
from symplyphysics import (
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    convert_to,
)

# Description
## The dimensionless form of the van der Waals equation of state features reduced quantities,
## which are simply the usual thermodynamic quantities divided by their value at the [critical
## point](TODO: add link). One notable property of the dimensionless equation of state is that
## it contains no substance-specific quantities, i.e. all van der Waals fluids will plot on the
## same reduced pressure-volume curve at the same reduced temperature.

# Law: (p* + 3 / (V*)**2) * (V* - 1/3) = 8/3 * T*
## p* = p/p_c - reduced pressure
## V* = V/V_c - reduced volume
## T* = T/T_c - reduced temperature

reduced_pressure = Symbol("reduced_pressure", dimensionless)
reduced_volume = Symbol("reduced_volume", dimensionless)
reduced_temperature = Symbol("reduced_temperature", dimensionless)

law = Eq(
    (reduced_pressure + 3 / reduced_volume**2) * (reduced_volume - Rational(1, 3)),
    Rational(8, 3) * reduced_temperature,
)

# TODO: derive from van der Waals equation of state


def print_law() -> str:
    return print_expression(law)


@validate_input(
    reduced_volume_=reduced_volume,
    reduced_temperature_=reduced_temperature,
)
@validate_output(reduced_pressure)
def calculate_reduced_pressure(
    reduced_volume_: float,
    reduced_temperature_: float,
) -> float:
    expr = solve(law, reduced_pressure)[0]
    result = expr.subs({
        reduced_volume: reduced_volume_,
        reduced_temperature: reduced_temperature_,
    })
    return float(convert_to(Quantity(result), S.One))
