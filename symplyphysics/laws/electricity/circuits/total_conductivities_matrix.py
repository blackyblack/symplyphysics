from sympy import Eq, solve, Matrix
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    SI,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

## Description
## The conductivity matrix is one of the ways to describe a microwave device. The Y-parameters of the device act as elements
## of this matrix. The matrix equation relates the input and output currents to the input and output voltages.

## Law is: [I1, I2] = Matrix([[Y11, Y12], [Y21, Y22]]) * [U1, U2], where
## U1 - input voltage,
## U2 - output voltage,
## Y11 - ratio of input current to input voltage in case of a short circuit at the output,
## Y12 - ratio of input current to output voltage in case of a short circuit at the output,
## Y21 - ratio of output current to input voltage in case of a short circuit at the output,
## Y22 - ratio of output current to output voltage in case of a short circuit at the output,
## I1 - input current,
## I2 - output current.

input_voltage = Symbol("input_voltage", units.voltage)
output_voltage = Symbol("output_voltage", units.voltage)
input_current = Symbol("input_current", units.current)
output_current = Symbol("output_current", units.current)
conductivity_input_input = Symbol("conductivity_input_input", units.conductance)
conductivity_input_output = Symbol("conductivity_input_output", units.conductance)
conductivity_output_input = Symbol("conductivity_output_input", units.conductance)
conductivity_output_output = Symbol("conductivity_output_output", units.conductance)

law = Eq(Matrix([input_current, output_current]), Matrix([[conductivity_input_input, conductivity_input_output], [conductivity_output_input, conductivity_output_output]]) * Matrix([input_voltage, output_voltage]))


def print_law() -> str:
    return print_expression(law)


@validate_input(input_current_=input_current, output_current_=output_current,)
@validate_output(units.voltage)
def calculate_voltages(input_current_: Quantity, output_current_: Quantity, conductivities_: tuple[tuple[Quantity, Quantity], tuple[Quantity, Quantity]]) -> tuple[Quantity, Quantity]:
    for conductivity in [x for y in conductivities_ for x in y]:
        assert_equivalent_dimension(conductivity, conductivity.dimension.name, "calculate_voltages", units.conductance)
    result_voltages = solve(law, [input_voltage, output_voltage], dict=True)[0]
    result_input_voltage = result_voltages[input_voltage]
    result_output_voltage = result_voltages[output_voltage]
    substitutions = {
        input_current: input_current_,
        output_current: output_current_,
        conductivity_input_input: conductivities_[0][0],
        conductivity_input_output: conductivities_[0][1],
        conductivity_output_input: conductivities_[1][0],
        conductivity_output_output: conductivities_[1][1],
    }
    result_input_voltage = result_input_voltage.subs(substitutions)
    result_output_voltage = result_output_voltage.subs(substitutions)
    return (Quantity(result_input_voltage), Quantity(result_output_voltage))
