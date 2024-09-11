"""
Oscillation period of inductor-capacitor network
================================================

An inductor-capacitor network, also called **LC circuit**, **resonator circuit**, or
**tuned circuit**, consists of an inductor and a capacitor connected together. This type
of circuit can act as an electrical resonator, storing energy oscillating at the circuit's
resonant frequency.
"""

from sympy import (Eq, Idx, solve, pi, sqrt, Derivative, simplify)
from symplyphysics import (units, Quantity, Symbol, Function, validate_input,
    validate_output, global_index)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.electricity import capacitance_from_charge_and_voltage as capacitance_definition
from symplyphysics.definitions import current_is_charge_derivative as charge_definition
from symplyphysics.definitions import harmonic_oscillator_is_second_derivative_equation as oscillator
from symplyphysics.laws.electricity import self_induced_electromotive_force_via_time_derivative_of_current as induction_voltage_definition
from symplyphysics.definitions import period_from_angular_frequency as period_definition
from symplyphysics.laws.electricity.circuits import sum_of_currents_through_junction_is_zero as kirchhoff_law
from symplyphysics.laws.electricity.circuits import sum_of_voltages_in_loop_is_zero as kirchhoff_law_2

period = Symbol("period", units.time)
"""
Natural period of oscillations.

Symbol:
    :code:`T`
"""

inductance = Symbol("inductance", units.inductance)
"""
Inductance of the inductor.

Symbol:
    :code:`L`
"""

capacitance = Symbol("capacitance", units.capacitance)
"""
Capacitance of the capacitor.

Symbol:
    :code:`C`
"""

law = Eq(period, 2 * pi * sqrt(inductance * capacitance))
r"""
:code:`T = 2 * pi * sqrt(L * C)`

Latex:
    .. math::
        T = 2 \pi \sqrt{L C}
"""

## Derive the same law from the capacitor, charge and self-induction voltage laws

## Let's assume we initially have capacitor charged to U0 voltage. In the zero time we connect this capacitor to inductor in a closed loop.
## So voltage on capacitor is always equals to voltage on inductor and the current through capacitor equals to current through inductor.

#NOTE: this proof is valid for capacitor and inductor in a closed loop without additional voltage source.
#      There are 2 more options: serial connection with external voltage source and parallel connection with external voltage source.
#      Additional proof can be added to show that oscillation period stays the same.

## 1. Prove that _capacitor_current(_time) = _inductor_current(_time)

_time = Symbol("_time", units.time)
_capacitor_current = Function("_capacitor_current", units.current)
_inductor_current = Function("_inductor_current", units.current)

_local_index_ = Idx("_local_index_", (1, 2))
_two_currents_law = kirchhoff_law.law.subs(global_index, _local_index_).doit()
# capacitor current is in, inductor current is out
_two_currents_applied = _two_currents_law.subs({
    kirchhoff_law.current[1]: _capacitor_current(_time),
    kirchhoff_law.current[2]: -1 * _inductor_current(_time)
})
_capacitor_current_applied = solve(_two_currents_applied, _capacitor_current(_time),
    dict=True)[0][_capacitor_current(_time)]
_capacitor_current_eq = Eq(_capacitor_current(_time), _capacitor_current_applied)

assert _capacitor_current_eq.lhs == _capacitor_current(_time)
assert _capacitor_current_eq.rhs == _inductor_current(_time)

## 2. Prove that _capacitor_voltage(_time) = _inductor_voltage(_time)

_capacitor_voltage = Function("_capacitor_voltage", units.voltage)
_inductor_voltage = Function("_inductor_voltage", units.voltage)

_two_voltages_law = kirchhoff_law_2.law.subs(global_index, _local_index_).doit()
# capacitor is voltage source, inductor is voltage consumer
_two_voltages_applied = _two_voltages_law.subs({
    kirchhoff_law_2.voltage[1]: -1 * _inductor_voltage(_time),
    kirchhoff_law_2.voltage[2]: _capacitor_voltage(_time)
})
_inductor_voltage_applied = solve(_two_voltages_applied, _inductor_voltage(_time),
    dict=True)[0][_inductor_voltage(_time)]

assert _inductor_voltage_applied == _capacitor_voltage(_time)

## 3. Prove that capacitor current derivative equals to capacitance * (second order derivative of voltage of capacitor)

## charge of capacitor is voltage of capacitor * capacitance
_capacitor_charge_law = capacitance_definition.definition.subs({
    capacitance_definition.capacitance: capacitance,
    capacitance_definition.charge: charge_definition.charge(_time),
    capacitance_definition.voltage: _capacitor_voltage(_time)
})
_capacitor_charge_applied = solve(_capacitor_charge_law, charge_definition.charge(_time),
    dict=True)[0][charge_definition.charge(_time)]

## I_c(t) = C * U_c'(t)
_capacitor_current_law = charge_definition.definition.subs(charge_definition.time, _time)
_capacitor_current_law = _capacitor_current_law.subs(charge_definition.charge(_time),
    _capacitor_charge_applied)
## I_c'(t) = C * U_c"(t)
_capacitor_current_law_derivative = Eq(Derivative(_capacitor_current_law.lhs, _time),
    Derivative(_capacitor_current_law.rhs, _time))

## 4. Prove that inductor voltage equals to -1 * capacitance * inductance * (second order derivative of voltage of capacitor)

## Inductor voltage is the self-inductance.
_inductor_voltage_law = induction_voltage_definition.law.subs(
    induction_voltage_definition.time, _time)
_inductor_voltage_law = _inductor_voltage_law.subs({
    induction_voltage_definition.inductance: inductance,
    induction_voltage_definition.electromotive_force(_time): _capacitor_voltage(_time),
    induction_voltage_definition.current(_time): charge_definition.current(_time)
})

_derived_law = [_inductor_voltage_law, _capacitor_current_law_derivative]

## U"(t) = - 1/LC * U(t)
_capacitor_voltage_solved = solve(_derived_law,
    (Derivative(charge_definition.current(_time)), _capacitor_voltage(_time)),
    dict=True)[0][_capacitor_voltage(_time)]
_voltage_diff_eq = Eq(_capacitor_voltage(_time), _capacitor_voltage_solved)

## 5. Solve differential equation and find period of the harmonic oscillator

## Expected solution for U"(t) = - 1/LC * U(t) is:
## A * e^(i * w * t) + B * e^(-i * w * t), where w = 1 / sqrt(LC)

_oscillator_eq = oscillator.definition.subs(oscillator.time, _time)
_oscillator_eq = _oscillator_eq.subs(oscillator.displacement(_time), _capacitor_voltage(_time))
angular_frequency_solved = simplify(
    solve([_oscillator_eq, _voltage_diff_eq], (oscillator.angular_frequency, _capacitor_voltage(_time)),
    dict=True)[0][oscillator.angular_frequency])

# 6. Derive period from frequency
_period_law = period_definition.law.subs(period_definition.angular_frequency,
    angular_frequency_solved)
_period_solved = solve(_period_law, period_definition.period, dict=True)[0][period_definition.period]
## Square roots fail to compare with each other. Raise both parts to power of 2 before checking for equality.
assert expr_equals(_period_solved**2, law.rhs**2)


@validate_input(inductance_=inductance, capacitance_=capacitance)
@validate_output(period)
def calculate_oscillation_period(inductance_: Quantity, capacitance_: Quantity) -> Quantity:
    result_period_expr = solve(law, period, dict=True)[0][period]
    result_expr = result_period_expr.subs({inductance: inductance_, capacitance: capacitance_})
    return Quantity(result_expr)
