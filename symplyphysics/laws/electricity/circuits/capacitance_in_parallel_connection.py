"""
Capacitance in parallel connection
==================================

The total capacitance of capacitors connected in parallel is the sum of the
capacitances of individual capacitors.

**Conditions:**

#. Components are connected in parallel.
"""

from sympy import Eq, Idx, solve, symbols as sympy_symbols
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    SumIndexed,
    global_index,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.symbols import clone_as_indexed
from symplyphysics.laws.electricity import (
    capacitance_from_charge_and_voltage as _capacitance_law,)

total_capacitance = symbols.capacitance
"""
Total :symbols:`capacitance`.
"""

capacitance = clone_as_indexed(symbols.capacitance)
"""
:symbols:`capacitance` of :math:`i`-th capacitor.
"""

law = Eq(total_capacitance, SumIndexed(capacitance[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""

# Derive from expression of capacitance via charge and voltage.

# Voltage is the same for all components connected in parallel.
# The proof is for two capacitors, but the same can be done for more capacitors.

_voltage = sympy_symbols("voltage")

_charge_expr = solve(_capacitance_law.definition,
    _capacitance_law.charge)[0].subs(_capacitance_law.voltage, _voltage)

_charge1 = _charge_expr.subs(_capacitance_law.capacitance, capacitance[1])
_charge2 = _charge_expr.subs(_capacitance_law.capacitance, capacitance[2])

# TODO Create law of additivity of electric charge
_total_charge_expr = _charge1 + _charge2

_total_capacitance_expr = _capacitance_law.definition.rhs.subs({
    _capacitance_law.charge: _total_charge_expr,
    _capacitance_law.voltage: _voltage,
})

_local_idx = Idx("i", (1, 2))
_capacitance_from_law = law.rhs.subs(global_index, _local_idx).doit()

assert expr_equals(_total_capacitance_expr, _capacitance_from_law)


@validate_input(capacitances_=capacitance)
@validate_output(total_capacitance)
def calculate_parallel_capacitance(capacitances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(capacitances_)))
    capacitances_law = law.subs(global_index, local_index)
    capacitances_law = capacitances_law.doit()
    solved = solve(capacitances_law, total_capacitance, dict=True)[0][total_capacitance]
    for i, v in enumerate(capacitances_):
        solved = solved.subs(capacitance[i + 1], v)
    return Quantity(solved)
