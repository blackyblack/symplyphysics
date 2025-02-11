"""
Impedances of T-type circuit of transmission line
=================================================

The T-type circuit consists of the first impedance connected in series, the third
impedance connected in parallel, and the second impedance connected in series. For a
transmission line, these impedances can be calculated for an equivalent replacement
circuit by knowing the transmission line parameters.

..
    TODO: find link
"""

from sympy import Eq, solve, Matrix, I, sinh, tanh, cosh, Symbol as SymSymbol, S
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity.circuits.transmission_lines import transmission_matrix_lossy_transmission_line as matrix_lossy_law
from symplyphysics.laws.electricity.circuits.transmission_lines import transmission_matrix_of_t_type_circuit as matrix_circuit_law

first_impedance = clone_as_symbol(symbols.electrical_impedance, subscript="1")
"""
First :symbols:`electrical_impedance`.
"""

second_impedance = clone_as_symbol(symbols.electrical_impedance, subscript="2")
"""
Second :symbols:`electrical_impedance`.
"""

third_impedance = clone_as_symbol(symbols.electrical_impedance, subscript="3")
"""
Third :symbols:`electrical_impedance`.
"""

surge_impedance = symbols.surge_impedance
"""
:symbols:`surge_impedance` of the transmission line.
"""

length = symbols.length
"""
:symbols:`length` of the transmission line.
"""

propagation_constant = SymbolNew("b", 1 / units.length)
"""
The **propagation constant** is the inverse of the signal :symbols:`wavelength`.
"""

loss_factor = SymbolNew("a", 1 / units.length)
"""
The **loss factor** shows how many times the transmitted signal weakens per unit
:symbols:`length` of the transmission line.
"""

_expression = (loss_factor + I * propagation_constant) * length

law = Eq(
    Matrix([first_impedance, second_impedance, third_impedance]),
    surge_impedance * Matrix([
        tanh(S.Half * _expression),
        tanh(S.Half * _expression),
        1 / sinh(_expression),
    ]))
"""
..
    NOTE: SymPy still evaluates matrix multiplication

:laws:symbol::

:laws:latex::
"""

# This law might be derived via "transmission_matrix_lossy_transmission_line" law and
# "transmission_matrix_of_t_type_circuit" law.

_matrix_lossy_law_applied = matrix_lossy_law.law.subs({
    matrix_lossy_law.surge_impedance: surge_impedance,
    matrix_lossy_law.length: length,
    matrix_lossy_law.propagation_constant: propagation_constant,
    matrix_lossy_law.loss_factor: loss_factor,
})
_matrix_derived_1 = solve(_matrix_lossy_law_applied, [
    matrix_lossy_law.voltage_voltage_parameter, matrix_lossy_law.voltage_current_parameter,
    matrix_lossy_law.current_voltage_parameter, matrix_lossy_law.current_current_parameter
],
    dict=True)[0]

_matrix_circuit_law_applied = matrix_circuit_law.law.subs({
    matrix_circuit_law.voltage_voltage_parameter:
    _matrix_derived_1[matrix_lossy_law.voltage_voltage_parameter],
    matrix_circuit_law.voltage_current_parameter:
    _matrix_derived_1[matrix_lossy_law.voltage_current_parameter],
    matrix_circuit_law.current_voltage_parameter:
    _matrix_derived_1[matrix_lossy_law.current_voltage_parameter],
    matrix_circuit_law.current_current_parameter:
    _matrix_derived_1[matrix_lossy_law.current_current_parameter],
})

# HACK: SymPy is unable to solve it automatically. Try with known solution and check for zero output
_matrix_circuit_law_solved = _matrix_circuit_law_applied.subs({
    matrix_circuit_law.first_impedance: law.rhs[0],
    matrix_circuit_law.second_impedance: law.rhs[1],
    matrix_circuit_law.third_impedance: law.rhs[2],
})

# Help SymPy by replacing all instances of `_expression` to single symbol
_expression_symbol = SymSymbol("x")
_matrix_circuit_law_solved = _matrix_circuit_law_solved.subs({
    _expression: 2 * _expression_symbol,
    _expression / 2: _expression_symbol,
    _expression.expand(): 2 * _expression_symbol,
    _expression.expand() / 2: _expression_symbol,
})
# From `tanh` definition
_matrix_circuit_law_solved = _matrix_circuit_law_solved.subs(
    tanh(_expression_symbol),
    sinh(_expression_symbol) / cosh(_expression_symbol))

# Check if solution results in all zeroes.
assert expr_equals(_matrix_circuit_law_solved.lhs[0], _matrix_circuit_law_solved.rhs[0])
assert expr_equals(_matrix_circuit_law_solved.lhs[1], _matrix_circuit_law_solved.rhs[1])
assert expr_equals(_matrix_circuit_law_solved.lhs[2], _matrix_circuit_law_solved.rhs[2])


@validate_input(characteristic_resistance_=surge_impedance,
    line_length_=length,
    constant_propagation_=propagation_constant,
    loss_factor_=loss_factor)
@validate_output(units.impedance)
def calculate_impedances(characteristic_resistance_: Quantity, line_length_: Quantity,
    constant_propagation_: Quantity, loss_factor_: Quantity) -> tuple[Quantity, Quantity, Quantity]:
    result = solve(law, [first_impedance, second_impedance, third_impedance], dict=True)[0]
    result_a = result[first_impedance]
    result_b = result[second_impedance]
    result_c = result[third_impedance]
    substitutions = {
        surge_impedance: characteristic_resistance_,
        length: line_length_,
        propagation_constant: constant_propagation_,
        loss_factor: loss_factor_,
    }
    result_a = Quantity(result_a.subs(substitutions))
    result_b = Quantity(result_b.subs(substitutions))
    result_c = Quantity(result_c.subs(substitutions))
    return (result_a, result_b, result_c)
