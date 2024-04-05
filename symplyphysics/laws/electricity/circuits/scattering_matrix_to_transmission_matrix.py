from sympy import Eq, solve, Matrix
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    dimensionless
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

## Description
##

parameter_voltage_to_voltage = Symbol("parameter_voltage_to_voltage", dimensionless)
parameter_impedance = Symbol("parameter_impedance", units.impedance)
parameter_conductance = Symbol("parameter_conductance", units.conductance)
parameter_current_to_current = Symbol("parameter_current_to_current", dimensionless)
characteristic_resistance = Symbol("characteristic_resistance", units.impedance)

reflection_coefficient = Symbol("reflection_coefficient", dimensionless)
reverse_transmission_ratio = Symbol("reverse_transmission_ratio", dimensionless)
transmission_ratio = Symbol("transmission_ratio", dimensionless)
output_reflection_coefficient = Symbol("output_reflection_coefficient", dimensionless)

expression_A = ((1 + reflection_coefficient) * (1 - output_reflection_coefficient) + reverse_transmission_ratio * transmission_ratio) / (2 * transmission_ratio)
expression_B = characteristic_resistance * ((1 + reflection_coefficient) * (1 + output_reflection_coefficient) - reverse_transmission_ratio * transmission_ratio) / (2 * transmission_ratio)
expression_C = (1 / characteristic_resistance) * ((1 - reflection_coefficient) * (1 - output_reflection_coefficient) - reverse_transmission_ratio * transmission_ratio) / (2 * transmission_ratio)
expression_D = ((1 - reflection_coefficient) * (1 + output_reflection_coefficient) + reverse_transmission_ratio * transmission_ratio) / (2 * transmission_ratio)

law = Eq(Matrix([[parameter_voltage_to_voltage, parameter_impedance], [parameter_conductance, parameter_current_to_current]]),
         Matrix([[expression_A, expression_B], [expression_C, expression_D]]))


def print_law() -> str:
    return print_expression(law)


@validate_input(characteristic_resistance_=characteristic_resistance,)
def calculate_transmission_matrix(characteristic_resistance_: Quantity, parameters_: tuple[tuple[float, float], tuple[float, float]]) -> tuple[tuple[Quantity, float], tuple[float, Quantity]]:
    result = solve(law, [parameter_voltage_to_voltage, parameter_impedance, parameter_conductance, parameter_current_to_current], dict=True)[0]
    result_A = result[parameter_voltage_to_voltage]
    result_B = result[parameter_impedance]
    result_C = result[parameter_conductance]
    result_D = result[parameter_current_to_current]
    substitutions = {
        characteristic_resistance: characteristic_resistance_,
        reflection_coefficient: parameters_[0][0],
        reverse_transmission_ratio: parameters_[0][1],
        transmission_ratio: parameters_[1][0],
        output_reflection_coefficient: parameters_[1][1],
    }
    result_A = result_A.subs(substitutions)
    result_B = result_B.subs(substitutions)
    result_C = result_C.subs(substitutions)
    result_D = result_D.subs(substitutions)
    assert_equivalent_dimension(result_B, 'result_B', "calculate_transmission_matrix", units.impedance)
    assert_equivalent_dimension(result_C, 'result_C', "calculate_transmission_matrix", units.conductance)
    return ((Quantity(result_A), Quantity(result_B)), (Quantity(result_C), Quantity(result_D)))
