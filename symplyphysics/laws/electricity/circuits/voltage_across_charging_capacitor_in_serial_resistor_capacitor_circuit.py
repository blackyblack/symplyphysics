r"""
Voltage across charging capacitor in serial resistor-capacitor circuit
======================================================================

A series circuit containing only a resistor, a capacitor, a switch and a constant DC source of voltage
is known as a charging circuit. The capacitor is initially uncharged while the switch is open. The switch
is closed at :math:`t = 0`. Afterwards, the capacitor begins to charge and the voltage across it rises
while the voltage across the resistor begins to drop. Note that the voltage across the capacitor never
reaches the voltage of the source.

**Conditions:**

#. The circuit is a DC one.
"""

from sympy import (Derivative, Eq, Idx, solve, exp, simplify)
from symplyphysics import (units, Quantity, Symbol, Function, validate_input,
    validate_output, global_index)
from symplyphysics.definitions import current_is_charge_derivative as charge_definition
from symplyphysics.definitions import capacitance_from_charge_and_voltage as capacitance_definition
from symplyphysics.laws.electricity import current_is_voltage_over_resistance as ohms_law
from symplyphysics.laws.electricity.circuits import sum_of_currents_through_junction_is_zero as kirchhoff_law
from symplyphysics.laws.electricity.circuits import sum_of_voltages_in_loop_is_zero as kirchhoff_law_2
from symplyphysics.laws.electricity.circuits import time_constant_of_resistor_capacitor_circuit as time_constant_law

capacitor_voltage = Symbol("capacitor_voltage", units.voltage)
"""
Voltage across the capacitor.

Symbol:
    :code:`V`
"""

source_voltage = Symbol("source_voltage", units.voltage)
r"""
Voltage of the source, which is the initial voltage across the resistor.

Symbol:
    :code:`V_0`

Latex:
    :math:`V_0`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

time_constant = Symbol("time_constant", units.time)
r"""
:doc:`Time constant <laws.electricity.circuits.time_constant_of_resistor_capacitor_circuit>` of the circuit.

Symbol:
    :code:`tau`

Latex:
    :math:`\tau`
"""

law = Eq(capacitor_voltage, source_voltage * (1 - exp(-time / time_constant)))
r"""
:code:`V = V_0 * (1 - exp(-1 * t / tau))`

Latex:
    .. math::
        V = V_0 \left( 1 - \exp \left( - \frac{t}{\tau} \right) \right)
"""

## Derive the same law from the Ohms and Kirchhoff laws

_resistance = time_constant_law.resistance
_capacitance = time_constant_law.capacitance

_capacitor_current = Function("_capacitor_current", units.current)
_resistor_current = Function("_resistor_current", units.current)
_resistor_voltage = Function("_resistor_voltage", units.voltage)

_local_index_ = Idx("_local_index_", (1, 2))
_two_currents_law = kirchhoff_law.law.subs(global_index, _local_index_).doit()
# capacitor current is in, resistor current is out
_two_currents_applied = _two_currents_law.subs({
    kirchhoff_law.current[1]: _capacitor_current(time),
    kirchhoff_law.current[2]: -1 * _resistor_current(time)
})
_capacitor_current_applied = solve(_two_currents_applied, _capacitor_current(time),
    dict=True)[0][_capacitor_current(time)]
_capacitor_current_eq = Eq(_capacitor_current(time), _capacitor_current_applied)

## 1. Prove that _capacitor_current(time) = _resistor_current(time)
assert _capacitor_current_eq.lhs == _capacitor_current(time)
assert _capacitor_current_eq.rhs == _resistor_current(time)

_local_index_ = Idx("_local_index_", (1, 3))
_three_voltages_law = kirchhoff_law_2.law.subs(global_index, _local_index_).doit()
# source_voltage is voltage source, capacitor and resistor are voltage consumers
_three_voltages_applied = _three_voltages_law.subs({
    kirchhoff_law_2.voltage[1]: -1 * capacitor_voltage,
    kirchhoff_law_2.voltage[2]: -1 * _resistor_voltage(time),
    kirchhoff_law_2.voltage[3]: source_voltage
})
_resistor_voltage_applied = solve(_three_voltages_applied, _resistor_voltage(time),
    dict=True)[0][_resistor_voltage(time)]
_resistor_voltage_eq = Eq(_resistor_voltage(time), _resistor_voltage_applied)

## 2. Prove that _resistor_voltage(time) = source_voltage - capacitor_voltage(time)
assert _resistor_voltage_eq.rhs == source_voltage - capacitor_voltage

# use _resistor_voltage as proven in _resistor_voltage_eq
# use charge_definition.current since it is same on resistor and capacitor as proven in _capacitor_current_eq
_resistor_ohm_eq = ohms_law.law.subs({
    ohms_law.voltage: source_voltage - capacitor_voltage,
    ohms_law.resistance: _resistance,
    ohms_law.current: charge_definition.current(time)
})
_capacitance_eq = capacitance_definition.definition.subs({
    capacitance_definition.capacitance: _capacitance,
    capacitance_definition.charge: charge_definition.charge(time),
    capacitance_definition.voltage: capacitor_voltage
})
_charge_eq = charge_definition.definition.subs(charge_definition.time, time)

_derived_law = [_resistor_ohm_eq, _capacitance_eq, _charge_eq]

_solved_charge_function = solve(_derived_law,
    (capacitor_voltage, charge_definition.current(time), charge_definition.charge(time)),
    dict=True)[0][charge_definition.charge(time)]
_charge_diff_eq = Eq(charge_definition.charge(time), _solved_charge_function)

## 3. Prove that charge(time) = capacitance * source_voltage - capacitance * resistance * Derivative(charge(time), time))
## Q(t) = U0 * C - R * C * dQ(t) / dt
_capacitor_charge_function = source_voltage * _capacitance - _resistance * _capacitance * Derivative(
    charge_definition.charge(time), time)
assert simplify(_charge_diff_eq.rhs - _capacitor_charge_function) == 0

## 4. Convert charge to capacitor voltage
_capacitor_voltage_solved = solve(_capacitance_eq, charge_definition.charge(time),
    dict=True)[0][charge_definition.charge(time)]
_voltage_diff_eq = _charge_diff_eq.subs(charge_definition.charge(time), _capacitor_voltage_solved)

## 5. Solve differential equation
# HACK: use known solution since sympy.dsolve() gives us another result

_rhs = law.rhs.subs(time_constant, time_constant_law.law.rhs)
_voltage_diff_solution = _voltage_diff_eq.subs(capacitor_voltage, _rhs)
assert simplify(_voltage_diff_solution.lhs - _voltage_diff_solution.rhs) == 0


@validate_input(initial_voltage_=source_voltage,
    time_constant_=time_constant,
    time_=time)
@validate_output(capacitor_voltage)
def calculate_capacitor_voltage(initial_voltage_: Quantity, time_constant_: Quantity,
    time_: Quantity) -> Quantity:
    capacitor_voltage_expr = solve(law, capacitor_voltage,
        dict=True)[0][capacitor_voltage]
    result_expr = capacitor_voltage_expr.subs({
        source_voltage: initial_voltage_,
        time_constant: time_constant_,
        time: time_
    })
    return Quantity(result_expr)
