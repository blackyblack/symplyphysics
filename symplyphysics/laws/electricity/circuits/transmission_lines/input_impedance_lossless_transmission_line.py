from sympy import Eq, solve, tan, I
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

## Description
## Knowing the length of the transmission line and its characteristic resistance, as well as the propagation constant and
## load resistance, it is possible to determine the input impedance of the transmission line.
## The propagation constant is the inverse of the wavelength.

## Law is: Zin = Z0 * (Zl + I * Z0 * tan(b * l)) / (Z0 + I * Zl * tan(b * l)), where
## Zin - input impedance of the transmission line,
## Zl - load resistance,
## Z0 - characteristic resistance of the transmission line,
## l - length of the transmission line,
## b - constant propagation of signal,
## I - imaginary unit.

input_impedance = Symbol("input_impedance", units.impedance)

characteristic_resistance = Symbol("characteristic_resistance", units.impedance)
load_resistance = Symbol("load_resistance", units.impedance)
constant_propagation = Symbol("constant_propagation", 1 / units.length)
line_length = Symbol("line_length", units.length)

law = Eq(
    input_impedance,
    characteristic_resistance *
    (load_resistance + I * characteristic_resistance * tan(constant_propagation * line_length)) /
    (characteristic_resistance + I * load_resistance * tan(constant_propagation * line_length)))


def print_law() -> str:
    return print_expression(law)


@validate_input(characteristic_resistance_=characteristic_resistance,
    load_resistance_=load_resistance,
    constant_propagation_=constant_propagation,
    line_length_=line_length)
@validate_output(input_impedance)
def calculate_input_impedance(characteristic_resistance_: Quantity, load_resistance_: Quantity,
    constant_propagation_: Quantity, line_length_: Quantity) -> Quantity:
    result_expr = solve(law, input_impedance, dict=True)[0][input_impedance]
    result_expr = result_expr.subs({
        characteristic_resistance: characteristic_resistance_,
        load_resistance: load_resistance_,
        constant_propagation: constant_propagation_,
        line_length: line_length_,
    })
    return Quantity(result_expr)
