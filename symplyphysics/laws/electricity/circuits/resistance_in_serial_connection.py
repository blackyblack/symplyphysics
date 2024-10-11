"""
Resistivity of serial resistors
===============================

The total resistance of the circuit whose components are connected in series is the sum
of the resistances of individual components.

**Conditions:**

#. Applies to direct current circuits.
"""

from sympy import Eq, Idx, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    SumIndexed,
    global_index,
    symbols,
)
from symplyphysics.core.symbols.symbols import clone_as_indexed
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import (
    current_is_voltage_over_resistance as _ohm_law,)

total_resistance = symbols.electrical_resistance
"""
Total resistance of the circuit.
"""

resistance = clone_as_indexed(symbols.electrical_resistance,
    display_symbol="R[i]",
    display_latex="R_i")
r"""
Resistance of the :math:`i`-th component.
"""

law = Eq(total_resistance, SumIndexed(resistance[global_index], global_index))
r"""
:laws:symbol::

:laws:latex::
"""

# Derive for two resistors.
# In a parallel connection, the same current flows through the components.

_voltage_expr = solve(_ohm_law.law, _ohm_law.voltage)[0]

_voltage1 = _voltage_expr.subs(_ohm_law.resistance, resistance[1])

_voltage2 = _voltage_expr.subs(_ohm_law.resistance, resistance[2])

# TODO: create law of voltage in serial connection and use it here
_total_voltage = _voltage1 + _voltage2

_total_resistance = solve(_ohm_law.law, _ohm_law.resistance)[0].subs(
    _ohm_law.voltage,
    _total_voltage,
)

_local_idx = Idx("local_idx", (1, 2))
_resistance_from_law = law.rhs.subs(global_index, _local_idx).doit()

assert expr_equals(_total_resistance, _resistance_from_law)


@validate_input(resistances_=resistance)
@validate_output(total_resistance)
def calculate_serial_resistance(resistances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(resistances_)))
    resistances_law = law.subs(global_index, local_index)
    resistances_law = resistances_law.doit()
    solved = solve(resistances_law, total_resistance, dict=True)[0][total_resistance]
    for i, v in enumerate(resistances_):
        solved = solved.subs(resistance[i + 1], v)
    return Quantity(solved)
