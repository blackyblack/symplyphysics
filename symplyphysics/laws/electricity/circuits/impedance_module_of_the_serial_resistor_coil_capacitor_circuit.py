"""
Impedance module of serial resistor-coil-capacitor circuit
==========================================================

Consider an electrical circuit consisting of a capacitor, coil, and resistor connected
in series. Then you can find the impedance module of such a circuit.

**Links:**

#. `Wikipedia, derivable from first formula <https://en.wikipedia.org/wiki/Electrical_impedance#Resistance_vs_reactance>`__.
"""

from sympy import (Eq, solve, sqrt, Idx, expand)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, global_index, symbols, clone_as_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity.circuits import capacitor_impedance_from_capacitive_reactance as capacitor_impedance_law
from symplyphysics.laws.electricity.circuits import coil_impedance_from_inductive_reactance as coil_impedance_law
from symplyphysics.laws.electricity.circuits import impedance_in_serial_connection as serial_law

circuit_impedance_module = Symbol("abs(Z)", units.impedance, display_latex="|Z|")
"""
Absolute value of the circuit's :symbols:`electrical_impedance`.
"""

resistor_resistance = clone_as_symbol(symbols.electrical_resistance, real=True)
"""
:symbols:`electrical_resistance` of the resistor.
"""

capacitor_reactance = clone_as_symbol(symbols.electrical_reactance, display_symbol="X_C", display_latex="X_\\text{C}", real=True)
"""
:symbols:`electrical_reactance` of the capacitor.
"""

coil_reactance = clone_as_symbol(symbols.electrical_reactance, display_symbol="X_L", display_latex="X_\\text{L}", real=True)
"""
:symbols:`electrical_reactance` of the coil.
"""

law = Eq(circuit_impedance_module,
    sqrt(resistor_resistance**2 + (coil_reactance - capacitor_reactance)**2))
"""
:laws:symbol::

:laws:latex::
"""

# This law might be derived via "serial_impedance" law, "capacitor_impedance_from_capacitive_reactance" law,
# "coil_impedance_from_inductive_reactance" law.

_impedance_law_applied_1 = capacitor_impedance_law.law.subs({
    capacitor_impedance_law.reactance: capacitor_reactance,
})
_capacitive_impedance_derived = solve(_impedance_law_applied_1,
    capacitor_impedance_law.impedance,
    dict=True)[0][capacitor_impedance_law.impedance]

_impedance_law_applied_2 = coil_impedance_law.law.subs({
    coil_impedance_law.reactance: coil_reactance,
})
_coil_impedance_derived = solve(_impedance_law_applied_2,
    coil_impedance_law.impedance,
    dict=True)[0][coil_impedance_law.impedance]

_local_index = Idx("index_local", (1, 3))
_serial_law_applied = serial_law.law.subs(global_index, _local_index)
_serial_law_applied = _serial_law_applied.doit()
_circuit_impedance_derived = solve(_serial_law_applied, serial_law.total_impedance,
    dict=True)[0][serial_law.total_impedance]
for i, v in enumerate(
    (resistor_resistance, _coil_impedance_derived, _capacitive_impedance_derived)):
    _circuit_impedance_derived = _circuit_impedance_derived.subs(serial_law.impedance[i + 1], v)

assert expr_equals(abs(_circuit_impedance_derived), expand(law.rhs))


@validate_input(resistance_resistor_=resistor_resistance,
    capacitive_reactance_=capacitor_reactance,
    coil_reactance_=coil_reactance)
@validate_output(circuit_impedance_module)
def calculate_circuit_impedance_module(resistance_resistor_: Quantity,
    capacitive_reactance_: Quantity, coil_reactance_: Quantity) -> Quantity:
    result_expr = solve(law, circuit_impedance_module, dict=True)[0][circuit_impedance_module]
    result_expr = result_expr.subs({
        resistor_resistance: resistance_resistor_,
        capacitor_reactance: capacitive_reactance_,
        coil_reactance: coil_reactance_
    })
    return Quantity(result_expr)
