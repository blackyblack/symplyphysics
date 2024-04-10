from sympy import Eq, solve, Matrix, S
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    dimensionless,
    convert_to
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

## Description
## The transmission parameters matrix is one of the ways to describe a microwave device. The ABCD-parameters of the device act as elements
## of this matrix. The matrix equation relates the input voltage and input current to the output voltage and output current.
## Knowing the impedance of the load connected in series to the transmission line,
## it is possible to calculate the parameters A, B, C, D of the transmission matrix of the load.

## Law is: Matrix([[A, B], [C, D]]) = Matrix([[1, Zl], [0, 1]]), where
## A - ratio of input voltage to output voltage at idle at the output,
## B - ratio of input voltage to output current in case of a short circuit at the output,
## C - ratio of input current to output voltage at idle at the output,
## D - ratio of input current to output current in case of a short circuit at the output,
## Zl - the impedance of the load connected in series to the transmission line.

parameter_voltage_to_voltage = Symbol("parameter_voltage_to_voltage", dimensionless)
parameter_impedance = Symbol("parameter_impedance", units.impedance)
parameter_conductance = Symbol("parameter_conductance", units.conductance)
parameter_current_to_current = Symbol("parameter_current_to_current", dimensionless)

load_impedance = Symbol("load_impedance", units.impedance)

law = Eq(Matrix([[parameter_voltage_to_voltage, parameter_impedance], [parameter_conductance, parameter_current_to_current]]),
         Matrix([[1, load_impedance], [Quantity(0 * units.siemens), 1]]))


def print_law() -> str:
    return print_expression(law)


@validate_input(load_impedance_=load_impedance)
def calculate_transmission_matrix(load_impedance_: Quantity) -> tuple[tuple[float, Quantity], tuple[Quantity, float]]:
    result = solve(law, [parameter_voltage_to_voltage, parameter_impedance, parameter_conductance, parameter_current_to_current], dict=True)[0]
    result_A = result[parameter_voltage_to_voltage]
    result_B = result[parameter_impedance]
    result_C = result[parameter_conductance]
    result_D = result[parameter_current_to_current]
    substitutions = {
        load_impedance: load_impedance_,
    }
    result_A = float(convert_to(Quantity(result_A.subs(substitutions)), S.One).evalf())
    result_B = Quantity(result_B.subs(substitutions))
    result_C = Quantity(result_C.subs(substitutions))
    result_D = float(convert_to(Quantity(result_D.subs(substitutions)), S.One).evalf())
    assert_equivalent_dimension(result_B, 'result_B', "calculate_transmission_matrix", units.impedance)
    assert_equivalent_dimension(result_C, 'result_C', "calculate_transmission_matrix", units.conductance)
    return ((result_A, result_B), (result_C, result_D))
