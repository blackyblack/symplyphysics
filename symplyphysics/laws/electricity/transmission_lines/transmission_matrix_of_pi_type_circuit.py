from sympy import Eq, solve, Matrix
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    dimensionless,
    convert_to_float
)
from symplyphysics.core.dimensions import assert_equivalent_dimension
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity.transmission_lines import transmission_matrix_for_a_series_load_in_line as series_law
from symplyphysics.laws.electricity.transmission_lines import transmission_matrix_for_a_parallel_load_in_line as parallel_law

## Description
## The transmission parameters matrix is one of the ways to describe a microwave device. The ABCD-parameters of the device act as elements
## of this matrix. The matrix equation relates the input voltage and input current to the output voltage and output current.
## The π-type circuit consists of the first impedance connected in parallel, the third impedance connected in series, and the
## second impedance connected in parallel.
## Knowing impedances, it is possible to calculate the parameters A, B, C, D of the transmission matrix of this line.

## Law is: Matrix([[A, B], [C, D]]) = Matrix([[1 + Z3 / Z2, Z3], [(1 / Z1) + (1 / Z2) + Z3 / (Z1 * Z2), 1 + Z3 / Z1]]), where
## A - ratio of input voltage to output voltage at idle at the output,
## B - ratio of input voltage to output current in case of a short circuit at the output,
## C - ratio of input current to output voltage at idle at the output,
## D - ratio of input current to output current in case of a short circuit at the output,
## Z1 - first impedance,
## Z2 - second impedance,
## Z3 - third impedance.

parameter_voltage_to_voltage = Symbol("parameter_voltage_to_voltage", dimensionless)
parameter_impedance = Symbol("parameter_impedance", units.impedance)
parameter_conductance = Symbol("parameter_conductance", units.conductance)
parameter_current_to_current = Symbol("parameter_current_to_current", dimensionless)

first_impedance = Symbol("first_impedance", units.impedance)
second_impedance = Symbol("second_impedance", units.impedance)
third_impedance = Symbol("third_impedance", units.impedance)


law = Eq(Matrix([[parameter_voltage_to_voltage, parameter_impedance], [parameter_conductance, parameter_current_to_current]]),
         Matrix([[1 + third_impedance / second_impedance, third_impedance], [(1 / first_impedance) + (1 / second_impedance) + third_impedance / (first_impedance * second_impedance), 1 + third_impedance / first_impedance]]))

# This law might be derived via "transmission_matrix_for_a_series_load_in_line" law and
# "transmission_matrix_for_a_parallel_load_in_line" law.

parallel_law_applied_1 = parallel_law.law.subs({
    parallel_law.load_impedance: first_impedance,
})
matrix_derived_1 = solve(parallel_law_applied_1, [parallel_law.parameter_voltage_to_voltage, parallel_law.parameter_impedance, parallel_law.parameter_conductance, parallel_law.parameter_current_to_current], dict=True)[0]
matrix_derived_1 = Matrix([[matrix_derived_1[parallel_law.parameter_voltage_to_voltage], matrix_derived_1[parallel_law.parameter_impedance]], [matrix_derived_1[parallel_law.parameter_conductance], matrix_derived_1[parallel_law.parameter_current_to_current]]])

series_law_applied = series_law.law.subs({
    series_law.load_impedance: third_impedance,
})
matrix_derived_3 = solve(series_law_applied, [series_law.parameter_voltage_to_voltage, series_law.parameter_impedance, series_law.parameter_conductance, series_law.parameter_current_to_current], dict=True)[0]
matrix_derived_3 = Matrix([[matrix_derived_3[series_law.parameter_voltage_to_voltage], matrix_derived_3[series_law.parameter_impedance]], [matrix_derived_3[series_law.parameter_conductance], matrix_derived_3[series_law.parameter_current_to_current]]])

parallel_law_applied_2 = parallel_law.law.subs({
    parallel_law.load_impedance: second_impedance,
})
matrix_derived_2 = solve(parallel_law_applied_2, [parallel_law.parameter_voltage_to_voltage, parallel_law.parameter_impedance, parallel_law.parameter_conductance, parallel_law.parameter_current_to_current], dict=True)[0]
matrix_derived_2 = Matrix([[matrix_derived_2[parallel_law.parameter_voltage_to_voltage], matrix_derived_2[parallel_law.parameter_impedance]], [matrix_derived_2[parallel_law.parameter_conductance], matrix_derived_2[parallel_law.parameter_current_to_current]]])

# When cascading elements, the final transfer matrix will be equal to the product of the matrices of each element.
matrix_derived = matrix_derived_1 * matrix_derived_3 * matrix_derived_2

# Check if derived ABCD-parameters are same as declared.
assert expr_equals(matrix_derived[0, 0], law.rhs[0, 0])
assert expr_equals(matrix_derived[0, 1], law.rhs[0, 1])
assert expr_equals(matrix_derived[1, 0], law.rhs[1, 0])
assert expr_equals(matrix_derived[1, 1], law.rhs[1, 1])


def print_law() -> str:
    return print_expression(law)


@validate_input(impedances_=units.impedance)
def calculate_transmission_matrix(impedances_: tuple[Quantity, Quantity, Quantity]) -> tuple[tuple[float, Quantity], tuple[Quantity, float]]:
    result = solve(law, [parameter_voltage_to_voltage, parameter_impedance, parameter_conductance, parameter_current_to_current], dict=True)[0]
    result_A = result[parameter_voltage_to_voltage]
    result_B = result[parameter_impedance]
    result_C = result[parameter_conductance]
    result_D = result[parameter_current_to_current]
    substitutions = {
        first_impedance: impedances_[0],
        second_impedance: impedances_[1],
        third_impedance: impedances_[2],
    }
    result_A = convert_to_float(Quantity(result_A.subs(substitutions)))
    result_B = Quantity(result_B.subs(substitutions))
    result_C = Quantity(result_C.subs(substitutions))
    result_D = convert_to_float(Quantity(result_D.subs(substitutions)))
    assert_equivalent_dimension(result_B, 'result_B', "calculate_transmission_matrix", units.impedance)
    assert_equivalent_dimension(result_C, 'result_C', "calculate_transmission_matrix", units.conductance)
    return ((result_A, result_B), (result_C, result_D))
