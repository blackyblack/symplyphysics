"""
Impedance in serial connection
==============================

The total impedance of a circuit whose components are connected in series is the sum
of the impedances of individual components.

**Conditions:**

#. Components are connected in series.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Electrical_impedance#Series_combination>`__.
"""

from sympy import Eq, Idx, solve
from symplyphysics import (Quantity, validate_input, validate_output, IndexedSum, global_index,
    symbols)
from symplyphysics.core.symbols.symbols import clone_as_indexed
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.electricity import current_is_voltage_over_impedance as _ohms_law

total_impedance = symbols.electrical_impedance
"""
Total :symbols:`electrical_impedance` of the circuit.
"""

impedance = clone_as_indexed(symbols.electrical_impedance)
"""
:symbols:`electrical_impedance` of the :math:`i`-th component.
"""

law = Eq(total_impedance, IndexedSum(impedance[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""

# Derive from Ohm's law and the additive property of voltage for a two-component system. For more
# components the derivation is similar to this one.

# NOTE: in a serial connection the current flowing through the components is the same.
_voltage_expr = solve(_ohms_law.law, _ohms_law.voltage)[0]

_first_voltage_expr = _voltage_expr.subs(_ohms_law.impedance, impedance[1])
_second_voltage_expr = _voltage_expr.subs(_ohms_law.impedance, impedance[2])

# NOTE: voltage is additive.
_total_voltage_original = _first_voltage_expr + _second_voltage_expr

# We replace the two components with a single component that yields the same voltage and current.
_total_voltage_replaced = _voltage_expr.subs(_ohms_law.impedance, total_impedance)

_total_impedance_derived = solve(
    Eq(_total_voltage_original, _total_voltage_replaced),
    total_impedance,
)[0]

_total_impedance_expected = law.rhs.subs(global_index, Idx("i", (1, 2)))

assert expr_equals(_total_impedance_derived, _total_impedance_expected)


@validate_input(impedances_=impedance)
@validate_output(total_impedance)
def calculate_serial_impedance(impedances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(impedances_)))
    impedances_law = law.subs(global_index, local_index)
    impedances_law = impedances_law.doit()
    solved = solve(impedances_law, total_impedance, dict=True)[0][total_impedance]
    for i, v in enumerate(impedances_):
        solved = solved.subs(impedance[i + 1], v)
    return Quantity(solved)
