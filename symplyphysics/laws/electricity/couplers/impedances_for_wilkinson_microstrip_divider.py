from sympy import Eq, solve, Matrix, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)

## Description
## The Wilkinson divider is a device designed to divide the power of a microwave signal into two output ports.
## Different sections of the divider consist of a microstrip line of different widths. There are four such sections in total and each has its own impedance.
## https://habrastorage.org/getpro/habr/upload_files/c24/031/52e/c2403152e2b320ab1c4c44f970dee1f2.gif

## Law is: Matrix([Z1, Z2, Z3, Z4]) = Matrix([Z0 * sqrt(k(1 + k^2)), Z0 * sqrt((1 + k^2) / k^3), Z0 * sqrt(k), Z0 / sqrt(k)]), where
## Z1 - first impedance,
## Z2 - second impedance,
## Z3 - third impedance,
## Z4 - fourth impedance,
## Z0 - resistance of the transmission line to which the divider is connected,
## k - ratio coefficient of the power at the outputs of the divider.

first_impedance = Symbol("first_impedance", units.impedance)
second_impedance = Symbol("second_impedance", units.impedance)
third_impedance = Symbol("third_impedance", units.impedance)
fourth_impedance = Symbol("fourth_impedance", units.impedance)

characteristic_resistance = Symbol("characteristic_resistance", units.impedance)
ratio_of_power = Symbol("ratio_of_power", dimensionless)

law = Eq(
    Matrix([first_impedance, second_impedance, third_impedance, fourth_impedance]),
    Matrix([
    characteristic_resistance * sqrt(ratio_of_power * (1 + ratio_of_power**2)),
    characteristic_resistance * sqrt(
    (1 + ratio_of_power**2) / ratio_of_power**3), characteristic_resistance * sqrt(ratio_of_power),
    characteristic_resistance / sqrt(ratio_of_power)
    ]))


def print_law() -> str:
    return print_expression(law)


@validate_input(characteristic_resistance_=characteristic_resistance,
    ratio_of_power_=ratio_of_power)
@validate_output(units.impedance)
def calculate_impedances(characteristic_resistance_: Quantity,
    ratio_of_power_: float) -> tuple[Quantity, Quantity, Quantity, Quantity]:
    result = solve(law, [first_impedance, second_impedance, third_impedance, fourth_impedance],
        dict=True)[0]
    result_z1 = result[first_impedance]
    result_z2 = result[second_impedance]
    result_z3 = result[third_impedance]
    result_z4 = result[fourth_impedance]
    substitutions = {
        characteristic_resistance: characteristic_resistance_,
        ratio_of_power: ratio_of_power_,
    }
    result_z1 = Quantity(result_z1.subs(substitutions))
    result_z2 = Quantity(result_z2.subs(substitutions))
    result_z3 = Quantity(result_z3.subs(substitutions))
    result_z4 = Quantity(result_z4.subs(substitutions))
    return (result_z1, result_z2, result_z3, result_z4)
