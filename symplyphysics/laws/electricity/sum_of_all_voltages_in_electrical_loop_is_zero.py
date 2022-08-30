from typing import List
from sympy import Idx, IndexedBase, Sum
from symplyphysics.quantity_decorator import assert_equivalent_dimension
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## sum(U) = 0
## Where U is voltage applied to an element of an electrical loop.
## This law grows from definition of work of the potential field done to move a unit of charge in this field.
## This work is independent on a trajectory and depends only on start and end point.
## So if we move charge in electrical field and then go back to start point, zero summary job is done by field during thise movements.
## If we imagine that it's not right, we would be able to generate energy just by moving charge in electric field, it breaks the energy conservation law.
## Also this rule is applicable to units of mass moving in gravitation field.

voltage = IndexedBase('voltage')
voltages_total = symbols('voltages_total')
i = symbols('i', cls=Idx)

law = Eq(Sum(voltage[i], (i, 1, voltages_total)), 0)

def print():
    return pretty(law, use_unicode=False)

@validate_input(voltage_in = units.voltage)
@validate_output(units.voltage)
def calculate_current(voltage_in: Quantity) -> Quantity:
    two_voltages_law = law.subs(voltages_total, 2).doit()
    solved = solve(two_voltages_law, voltage[2], dict=True)[0][voltage[2]]
    result_expr = solved.subs(voltage[1], voltage_in)
    return expr_to_quantity(result_expr, 'voltage_out')

@validate_output(units.voltage)
def calculate_current_from_array(voltages: List[Quantity]) -> Quantity:
    for idx, c in enumerate(voltages):
        assert_equivalent_dimension(c, "validate_input", f"voltages[{idx}]", "calculate_voltage_from_array", units.voltage)
    voltages_law = law.subs(voltages_total, len(voltages) + 1).doit()
    unknown_voltage = voltage[len(voltages) + 1]
    solved = solve(voltages_law, unknown_voltage, dict=True)[0][unknown_voltage]
    for idx, c in enumerate(voltages):
        solved = solved.subs(voltage[idx + 1], c)
    return expr_to_quantity(solved, 'voltage_out')
