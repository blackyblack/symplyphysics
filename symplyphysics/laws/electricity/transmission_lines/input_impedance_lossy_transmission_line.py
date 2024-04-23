from sympy import Eq, solve, sinh, cosh, I
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

## Description
## Knowing the length of the transmission line, the loss factor of the transmission line and its characteristic resistance,
## as well as the propagation constant and load resistance, it is possible to determine the input impedance of the transmission line.
## The propagation constant is the inverse of the wavelength.
## The loss factor shows how many times the transmitted signal weakens per unit length of the transmission line.

## Law is: Zin = (cosh((a + I * b) * l) * Zl + Z0 * sinh((a + I * b) * l)) / (Zl * sinh((a + I * b) * l) / Z0 + cosh((a + I * b) * l)), where
## Zin - input impedance of the transmission line,
## Zl - load resistance,
## Z0 - characteristic resistance of the transmission line,
## l - length of the transmission line,
## b - constant propagation of signal,
## a - coefficient of signal loss in the transmission line,
## I - imaginary unit,
## sinh - hyperbolic sine,
## cosh - hyperbolic cosine.

input_impedance = Symbol("input_impedance", units.impedance)

characteristic_resistance = Symbol("characteristic_resistance", units.impedance)
load_resistance = Symbol("load_resistance", units.impedance)
constant_propagation = Symbol("constant_propagation", 1 / units.length)
line_length = Symbol("line_length", units.length)
loss_factor = Symbol("loss_factor", 1 / units.length)

expression = (loss_factor + I * constant_propagation) * line_length
law = Eq(input_impedance,
    (cosh(expression) * load_resistance + characteristic_resistance * sinh(expression)) /
    (load_resistance * sinh(expression) / characteristic_resistance + cosh(expression)))


def print_law() -> str:
    return print_expression(law)


@validate_input(characteristic_resistance_=characteristic_resistance,
    load_resistance_=load_resistance,
    constant_propagation_=constant_propagation,
    line_length_=line_length,
    loss_factor_=loss_factor)
@validate_output(input_impedance)
def calculate_input_impedance(characteristic_resistance_: Quantity, load_resistance_: Quantity,
    constant_propagation_: Quantity, line_length_: Quantity, loss_factor_: Quantity) -> Quantity:
    result_expr = solve(law, input_impedance, dict=True)[0][input_impedance]
    result_expr = result_expr.subs({
        characteristic_resistance: characteristic_resistance_,
        load_resistance: load_resistance_,
        constant_propagation: constant_propagation_,
        line_length: line_length_,
        loss_factor: loss_factor_,
    })
    return Quantity(result_expr)
