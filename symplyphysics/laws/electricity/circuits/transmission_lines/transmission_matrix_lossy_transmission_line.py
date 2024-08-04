from sympy import Eq, solve, Matrix, S, I, sinh, cosh
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    dimensionless,
    convert_to,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

## Description
## The transmission parameters matrix is one of the ways to describe a microwave device. The ABCD-parameters of the device act as elements
## of this matrix. The matrix equation relates the input voltage and input current to the output voltage and output current.
## Knowing the length and the loss factor of the transmission line, as well as the characteristic resistance of the line and the constant propagation of signal,
## it is possible to calculate the parameters A, B, C, D of the transmission matrix of this line.
## The propagation constant is the inverse of the wavelength.
## The loss factor shows how many times the transmitted signal weakens per unit length of the transmission line.

## Law is: Matrix([[A, B], [C, D]]) = Matrix([[cosh((a + I * b) * l), Z0 * sinh((a + I * b) * l)], [(1 / Z0) * sinh((a + I * b) * l), cosh((a + I * b) * l)]]), where
## A - ratio of input voltage to output voltage at idle at the output,
## B - ratio of input voltage to output current in case of a short circuit at the output,
## C - ratio of input current to output voltage at idle at the output,
## D - ratio of input current to output current in case of a short circuit at the output,
## Z0 - characteristic resistance of the transmission line,
## l - length of the transmission line,
## b - constant propagation of signal,
## a - coefficient of signal loss in the transmission line,
## I - imaginary unit,
## sinh - hyperbolic sine,
## cosh - hyperbolic cosine.

parameter_voltage_to_voltage = Symbol("parameter_voltage_to_voltage", dimensionless)
parameter_impedance = Symbol("parameter_impedance", units.impedance)
parameter_conductance = Symbol("parameter_conductance", units.conductance)
parameter_current_to_current = Symbol("parameter_current_to_current", dimensionless)

characteristic_resistance = Symbol("characteristic_resistance", units.impedance)
line_length = Symbol("line_length", units.length)
constant_propagation = Symbol("constant_propagation", 1 / units.length)
loss_factor = Symbol("loss_factor", 1 / units.length)

expression = (loss_factor + I * constant_propagation) * line_length
law = Eq(
    Matrix([[parameter_voltage_to_voltage, parameter_impedance],
    [parameter_conductance, parameter_current_to_current]]),
    Matrix([[cosh(expression), characteristic_resistance * sinh(expression)],
    [(1 / characteristic_resistance) * sinh(expression),
    cosh(expression)]]))


def print_law() -> str:
    return print_expression(law)


@validate_input(characteristic_resistance_=characteristic_resistance,
    line_length_=line_length,
    constant_propagation_=constant_propagation,
    loss_factor_=loss_factor)
def calculate_transmission_matrix(
        characteristic_resistance_: Quantity, line_length_: Quantity,
        constant_propagation_: Quantity,
        loss_factor_: Quantity) -> tuple[tuple[float, Quantity], tuple[Quantity, float]]:
    result = solve(law, [
        parameter_voltage_to_voltage, parameter_impedance, parameter_conductance,
        parameter_current_to_current
    ],
        dict=True)[0]
    result_a = result[parameter_voltage_to_voltage]
    result_b = result[parameter_impedance]
    result_c = result[parameter_conductance]
    result_d = result[parameter_current_to_current]
    substitutions = {
        characteristic_resistance: characteristic_resistance_,
        line_length: line_length_,
        constant_propagation: constant_propagation_,
        loss_factor: loss_factor_,
    }
    result_a = convert_to(result_a.subs(substitutions), S.One)
    result_b = Quantity(result_b.subs(substitutions))
    result_c = Quantity(result_c.subs(substitutions))
    result_d = convert_to(result_d.subs(substitutions), S.One)
    assert_equivalent_dimension(result_b, 'result_b', "calculate_transmission_matrix",
        units.impedance)
    assert_equivalent_dimension(result_c, 'result_c', "calculate_transmission_matrix",
        units.conductance)
    return ((result_a, result_b), (result_c, result_d))
