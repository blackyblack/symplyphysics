from sympy import Eq, solve, Matrix, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless
)

## Description
## The impedance matrix is one of the ways to describe a microwave device. The S-parameters of the device act as elements
## of this matrix. The matrix equation relates the waves reflected from the input and output to the waves incident on the
## input and output.
## a1, a2 - these coefficients are equal to the voltage incident on the port divided by the wave resistance of the transmission line.
## b1, b2 - these coefficients are equal to the voltage reflected from the port divided by the wave resistance of the transmission line.

## Law is: [b1, b2] = Matrix([[S11, S12], [S21, S22]]) * [a1, a2], where
## b1 - reflected wave coefficient per input port,
## b2 - reflected wave coefficient per output port,
## S11 - voltage reflection coefficient,
## S12 - reverse voltage transmission ratio,
## S21 - voltage transmission ratio,
## S22 - output voltage reflection coefficient,
## a1 - incident wave coefficient per input port,
## a2 - incident wave coefficient per output port.

input_reflected_wave = Symbol("input_reflected_wave", sqrt(units.power))
output_reflected_wave = Symbol("output_reflected_wave", sqrt(units.power))
input_wave = Symbol("input_wave", sqrt(units.power))
output_wave = Symbol("output_wave", sqrt(units.power))
reflection_coefficient = Symbol("reflection_coefficient", dimensionless)
reverse_transmission_ratio = Symbol("reverse_transmission_ratio", dimensionless)
transmission_ratio = Symbol("transmission_ratio", dimensionless)
output_reflection_coefficient = Symbol("output_reflection_coefficient", dimensionless)

law = Eq(Matrix([input_reflected_wave, output_reflected_wave]), Matrix([[reflection_coefficient, reverse_transmission_ratio], [transmission_ratio, output_reflection_coefficient]]) * Matrix([input_wave, output_wave]))


def print_law() -> str:
    return print_expression(law)


@validate_input(input_reflected_wave_=input_reflected_wave, output_reflected_wave_=output_reflected_wave)
@validate_output(sqrt(units.power))
def calculate_waves(input_reflected_wave_: Quantity, output_reflected_wave_: Quantity, parameters_: tuple[tuple[float, float], tuple[float, float]]) -> tuple[Quantity, Quantity]:
    result = solve(law, [input_wave, output_wave], dict=True)[0]
    result_input_wave = result[input_wave]
    result_output_wave = result[output_wave]
    substitutions = {
        input_reflected_wave: input_reflected_wave_,
        output_reflected_wave: output_reflected_wave_,
        reflection_coefficient: parameters_[0][0],
        reverse_transmission_ratio: parameters_[0][1],
        transmission_ratio: parameters_[1][0],
        output_reflection_coefficient: parameters_[1][1],
    }
    result_input_wave = result_input_wave.subs(substitutions)
    result_output_wave = result_output_wave.subs(substitutions)
    return (Quantity(result_input_wave), Quantity(result_output_wave))
