from sympy import Eq, solve, Matrix, I, sin, cos
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    dimensionless,
    convert_to_float,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

## Description
## The transmission parameters matrix is one of the ways to describe a microwave device. The ABCD-parameters of the device act as elements
## of this matrix. The matrix equation relates the input voltage and input current to the output voltage and output current.
## Knowing the length of the transmission line, as well as the characteristic resistance of the line and the constant propagation of signal,
## it is possible to calculate the parameters A, B, C, D of the transmission matrix of this line.
## The propagation constant is the inverse of the wavelength.

## Law is: Matrix([[A, B], [C, D]]) = Matrix([[cos(b * l), I * Z0 * sin(b * l)], [I * (1 / Z0) * sin(b * l), cos(b * l)]]), where
## A - ratio of input voltage to output voltage at idle at the output,
## B - ratio of input voltage to output current in case of a short circuit at the output,
## C - ratio of input current to output voltage at idle at the output,
## D - ratio of input current to output current in case of a short circuit at the output,
## Z0 - characteristic resistance of the transmission line,
## l - length of the transmission line,
## b - constant propagation of signal,
## I - imaginary unit.

parameter_voltage_to_voltage = Symbol("parameter_voltage_to_voltage", dimensionless)
parameter_impedance = Symbol("parameter_impedance", units.impedance)
parameter_conductance = Symbol("parameter_conductance", units.conductance)
parameter_current_to_current = Symbol("parameter_current_to_current", dimensionless)

characteristic_resistance = Symbol("characteristic_resistance", units.impedance)
line_length = Symbol("line_length", units.length)
constant_propagation = Symbol("constant_propagation", 1 / units.length)


law = Eq(Matrix([[parameter_voltage_to_voltage, parameter_impedance], [parameter_conductance, parameter_current_to_current]]),
         Matrix([[cos(constant_propagation * line_length), I * characteristic_resistance * sin(constant_propagation * line_length)], [I * (1 / characteristic_resistance) * sin(constant_propagation * line_length), cos(constant_propagation * line_length)]]))


def print_law() -> str:
    return print_expression(law)


@validate_input(characteristic_resistance_=characteristic_resistance, line_length_=line_length, constant_propagation_=constant_propagation)
def calculate_transmission_matrix(characteristic_resistance_: Quantity, line_length_: Quantity, constant_propagation_: Quantity) -> tuple[tuple[float, Quantity], tuple[Quantity, float]]:
    result = solve(law, [parameter_voltage_to_voltage, parameter_impedance, parameter_conductance, parameter_current_to_current], dict=True)[0]
    result_A = result[parameter_voltage_to_voltage]
    result_B = result[parameter_impedance]
    result_C = result[parameter_conductance]
    result_D = result[parameter_current_to_current]
    substitutions = {
        characteristic_resistance: characteristic_resistance_,
        line_length: line_length_,
        constant_propagation: constant_propagation_,
    }
    result_A = convert_to_float(Quantity(result_A.subs(substitutions)))
    result_B = Quantity(result_B.subs(substitutions))
    result_C = Quantity(result_C.subs(substitutions))
    result_D = convert_to_float(Quantity(result_D.subs(substitutions)))
    assert_equivalent_dimension(result_B, 'result_B', "calculate_transmission_matrix", units.impedance)
    assert_equivalent_dimension(result_C, 'result_C', "calculate_transmission_matrix", units.conductance)
    return ((result_A, result_B), (result_C, result_D))
