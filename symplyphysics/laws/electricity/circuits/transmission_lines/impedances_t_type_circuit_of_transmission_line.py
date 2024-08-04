from sympy import Eq, solve, Matrix, I, sinh, tanh, cosh, symbols
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity.circuits.transmission_lines import transmission_matrix_lossy_transmission_line as matrix_lossy_law
from symplyphysics.laws.electricity.circuits.transmission_lines import transmission_matrix_of_t_type_circuit as matrix_circuit_law

## Description
## The T-type circuit consists of the first impedance connected in series, the third impedance connected in parallel, and the
## second impedance connected in series.
## For a transmission line, these impedances can be calculated for an equivalent replacement circuit by knowing the transmission line parameters.
## The propagation constant is the inverse of the wavelength.
## The loss factor shows how many times the transmitted signal weakens per unit length of the transmission line.

## Law is: Matrix([Z1, Z2, Z3]) = Matrix([Z0 * tanh((a + I * b) * l / 2), Z0 * tanh((a + I * b) * l / 2), Z0 / sinh((a + I * b) * l)]), where
## Z1 - first impedance,
## Z2 - second impedance,
## Z3 - third impedance,
## Z0 - characteristic resistance of the transmission line,
## l - length of the transmission line,
## b - constant propagation of signal,
## a - coefficient of signal loss in the transmission line,
## I - imaginary unit,
## sinh - hyperbolic sine,
## tanh - hyperbolic tangent.

first_impedance = Symbol("first_impedance", units.impedance)
second_impedance = Symbol("second_impedance", units.impedance)
third_impedance = Symbol("third_impedance", units.impedance)

characteristic_resistance = Symbol("characteristic_resistance", units.impedance)
line_length = Symbol("line_length", units.length)
constant_propagation = Symbol("constant_propagation", 1 / units.length)
loss_factor = Symbol("loss_factor", 1 / units.length)

expression = (loss_factor + I * constant_propagation) * line_length
law = Eq(
    Matrix([first_impedance, second_impedance, third_impedance]),
    Matrix([
    characteristic_resistance * tanh(expression / 2),
    characteristic_resistance * tanh(expression / 2), characteristic_resistance / sinh(expression)
    ]))

# This law might be derived via "transmission_matrix_lossy_transmission_line" law and
# "transmission_matrix_of_t_type_circuit" law.

matrix_lossy_law_applied = matrix_lossy_law.law.subs({
    matrix_lossy_law.characteristic_resistance: characteristic_resistance,
    matrix_lossy_law.line_length: line_length,
    matrix_lossy_law.constant_propagation: constant_propagation,
    matrix_lossy_law.loss_factor: loss_factor,
})
matrix_derived_1 = solve(matrix_lossy_law_applied, [
    matrix_lossy_law.parameter_voltage_to_voltage, matrix_lossy_law.parameter_impedance,
    matrix_lossy_law.parameter_conductance, matrix_lossy_law.parameter_current_to_current
],
    dict=True)[0]

matrix_circuit_law_applied = matrix_circuit_law.law.subs({
    matrix_circuit_law.parameter_voltage_to_voltage:
    matrix_derived_1[matrix_lossy_law.parameter_voltage_to_voltage],
    matrix_circuit_law.parameter_impedance:
    matrix_derived_1[matrix_lossy_law.parameter_impedance],
    matrix_circuit_law.parameter_conductance:
    matrix_derived_1[matrix_lossy_law.parameter_conductance],
    matrix_circuit_law.parameter_current_to_current:
    matrix_derived_1[matrix_lossy_law.parameter_current_to_current],
})

# HACK: SymPy is unable to solve it automatically. Try with known solution and check for zero output
matrix_circuit_law_solved = matrix_circuit_law_applied.subs({
    matrix_circuit_law.first_impedance: law.rhs[0],
    matrix_circuit_law.second_impedance: law.rhs[1],
    matrix_circuit_law.third_impedance: law.rhs[2],
})

# Help SymPy by replacing all instances of `expression` to single symbol
expression_symbol = symbols("expression_symbol")
matrix_circuit_law_solved = matrix_circuit_law_solved.subs({
    expression: 2 * expression_symbol,
    expression / 2: expression_symbol,
    expression.expand(): 2 * expression_symbol,
    expression.expand() / 2: expression_symbol,
})
# From `tanh` definition
matrix_circuit_law_solved = matrix_circuit_law_solved.subs(
    tanh(expression_symbol),
    sinh(expression_symbol) / cosh(expression_symbol))

# Check if solution results in all zeroes.
assert expr_equals(matrix_circuit_law_solved.lhs[0], matrix_circuit_law_solved.rhs[0])
assert expr_equals(matrix_circuit_law_solved.lhs[1], matrix_circuit_law_solved.rhs[1])
assert expr_equals(matrix_circuit_law_solved.lhs[2], matrix_circuit_law_solved.rhs[2])


def print_law() -> str:
    return print_expression(law)


@validate_input(characteristic_resistance_=characteristic_resistance,
    line_length_=line_length,
    constant_propagation_=constant_propagation,
    loss_factor_=loss_factor)
@validate_output(units.impedance)
def calculate_impedances(characteristic_resistance_: Quantity, line_length_: Quantity,
    constant_propagation_: Quantity, loss_factor_: Quantity) -> tuple[Quantity, Quantity, Quantity]:
    result = solve(law, [first_impedance, second_impedance, third_impedance], dict=True)[0]
    result_a = result[first_impedance]
    result_b = result[second_impedance]
    result_c = result[third_impedance]
    substitutions = {
        characteristic_resistance: characteristic_resistance_,
        line_length: line_length_,
        constant_propagation: constant_propagation_,
        loss_factor: loss_factor_,
    }
    result_a = Quantity(result_a.subs(substitutions))
    result_b = Quantity(result_b.subs(substitutions))
    result_c = Quantity(result_c.subs(substitutions))
    return (result_a, result_b, result_c)
