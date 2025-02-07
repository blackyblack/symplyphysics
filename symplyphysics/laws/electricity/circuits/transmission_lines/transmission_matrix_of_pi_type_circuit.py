"""
Transmission matrix of π-type matrix
====================================

The π-type circuit consists of the first impedance connected in parallel, the third
impedance connected in series, and the second impedance connected in parallel. Knowing
the impedances, it is possible to calculate the parameters :math:`A, B, C, D` of the
transmission matrix of this line.

**Notes:**

#. See :ref:`Transmission matrix`.
#. Scheme of the circuit:

.. image:: https://en.wikipedia.org/wiki/%CE%A0_pad#/media/File:Attenuator,_Pi-section.svg
    :width: 400px
    :align: center

..
    TODO: find link
"""

from sympy import Eq, solve, Matrix
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    dimensionless,
    convert_to_float,
    symbols,
    clone_as_symbol,
    SymbolNew
)
from symplyphysics.core.dimensions import assert_equivalent_dimension
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity.circuits.transmission_lines import transmission_matrix_for_a_series_load_in_line as series_law
from symplyphysics.laws.electricity.circuits.transmission_lines import transmission_matrix_for_a_parallel_load_in_line as parallel_law

voltage_voltage_parameter = SymbolNew("A", dimensionless)
"""
Ratio of input :symbols:`voltage` to output :symbols:`voltage` at idle at the output.
"""

voltage_current_parameter = SymbolNew("B", units.impedance)
"""
Ratio of input :symbols:`voltage` to output :symbols:`current` in case of a short circuit at the output.
"""

current_voltage_parameter = SymbolNew("C", units.conductance)
"""
Ratio of input :symbols:`current` to output :symbols:`voltage` at idle at the output.
"""

current_current_parameter = SymbolNew("D", dimensionless)
"""
Ratio of input :symbols:`current` to output :symbols:`current` in case of a short circuit at the output.
"""

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

law = Eq(
    Matrix([[voltage_voltage_parameter, voltage_current_parameter], [current_voltage_parameter, current_current_parameter]]),
    Matrix([[1 + third_impedance / second_impedance, third_impedance],
    [(1 / first_impedance) + (1 / second_impedance) + third_impedance / (first_impedance * second_impedance), 1 + third_impedance / first_impedance]]))
"""
..
    NOTE: Code printing in upcoming PR.
"""

# This law might be derived via "transmission_matrix_for_a_series_load_in_line" law and
# "transmission_matrix_for_a_parallel_load_in_line" law.

_parallel_law_applied_1 = parallel_law.law.subs({
    parallel_law.load_impedance: first_impedance,
})
_matrix_derived_1 = solve(_parallel_law_applied_1, [
    parallel_law.voltage_voltage_parameter, parallel_law.voltage_current_parameter,
    parallel_law.current_voltage_parameter, parallel_law.current_current_parameter
],
    dict=True)[0]
_matrix_derived_1 = Matrix([[
    _matrix_derived_1[parallel_law.voltage_voltage_parameter],
    _matrix_derived_1[parallel_law.voltage_current_parameter]
],
    [
    _matrix_derived_1[parallel_law.current_voltage_parameter],
    _matrix_derived_1[parallel_law.current_current_parameter]
    ]])

_series_law_applied = series_law.law.subs({
    series_law.load_impedance: third_impedance,
})
_matrix_derived_3 = solve(_series_law_applied, [
    series_law.voltage_voltage_parameter, series_law.voltage_current_parameter,
    series_law.current_voltage_parameter, series_law.current_current_parameter
],
    dict=True)[0]
_matrix_derived_3 = Matrix([[
    _matrix_derived_3[series_law.voltage_voltage_parameter],
    _matrix_derived_3[series_law.voltage_current_parameter]
],
    [
    _matrix_derived_3[series_law.current_voltage_parameter],
    _matrix_derived_3[series_law.current_current_parameter]
    ]])

_parallel_law_applied_2 = parallel_law.law.subs({
    parallel_law.load_impedance: second_impedance,
})
_matrix_derived_2 = solve(_parallel_law_applied_2, [
    parallel_law.voltage_voltage_parameter, parallel_law.voltage_current_parameter,
    parallel_law.current_voltage_parameter, parallel_law.current_current_parameter
],
    dict=True)[0]
_matrix_derived_2 = Matrix([[
    _matrix_derived_2[parallel_law.voltage_voltage_parameter],
    _matrix_derived_2[parallel_law.voltage_current_parameter]
],
    [
    _matrix_derived_2[parallel_law.current_voltage_parameter],
    _matrix_derived_2[parallel_law.current_current_parameter]
    ]])

# When cascading elements, the final transfer matrix will be equal to the product of the matrices of each element.
_matrix_derived = _matrix_derived_1 * _matrix_derived_3 * _matrix_derived_2

# Check if derived ABCD-parameters are same as declared.
assert expr_equals(_matrix_derived[0, 0], law.rhs[0, 0])
assert expr_equals(_matrix_derived[0, 1], law.rhs[0, 1])
assert expr_equals(_matrix_derived[1, 0], law.rhs[1, 0])
assert expr_equals(_matrix_derived[1, 1], law.rhs[1, 1])


@validate_input(impedances_=units.impedance)
def calculate_transmission_matrix(
    impedances_: tuple[Quantity, Quantity, Quantity]
) -> tuple[tuple[float, Quantity], tuple[Quantity, float]]:
    result = solve(law, [
        voltage_voltage_parameter, voltage_current_parameter, current_voltage_parameter,
        current_current_parameter
    ],
        dict=True)[0]
    result_a = result[voltage_voltage_parameter]
    result_b = result[voltage_current_parameter]
    result_c = result[current_voltage_parameter]
    result_d = result[current_current_parameter]
    substitutions = {
        first_impedance: impedances_[0],
        second_impedance: impedances_[1],
        third_impedance: impedances_[2],
    }
    result_a = convert_to_float(Quantity(result_a.subs(substitutions)))
    result_b = Quantity(result_b.subs(substitutions))
    result_c = Quantity(result_c.subs(substitutions))
    result_d = convert_to_float(Quantity(result_d.subs(substitutions)))
    assert_equivalent_dimension(result_b, 'result_b', "calculate_transmission_matrix",
        units.impedance)
    assert_equivalent_dimension(result_c, 'result_c', "calculate_transmission_matrix",
        units.conductance)
    return ((result_a, result_b), (result_c, result_d))
