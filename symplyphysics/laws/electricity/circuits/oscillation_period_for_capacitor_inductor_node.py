from sympy import (Eq, Idx, solve, pi, sqrt, Derivative, simplify)
from symplyphysics import (units, Quantity, Symbol, Function, print_expression, validate_input,
    validate_output, global_index)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.definitions import capacitance_from_charge_and_voltage as capacitance_definition
from symplyphysics.definitions import current_is_charge_derivative as charge_definition
from symplyphysics.definitions import harmonic_oscillator_is_second_derivative_equation as oscillator
from symplyphysics.laws.electricity import self_induction_voltage_from_current_derivative as induction_voltage_definition
from symplyphysics.definitions import period_from_angular_frequency as period_definition
from symplyphysics.laws.electricity.circuits import sum_of_currents_through_junction_is_zero as kirchhoff_law
from symplyphysics.laws.electricity.circuits import sum_of_all_voltages_in_loop_is_zero as kirchhoff_law_2

# Description
## LC-oscillator is the circuit of inductor and capacitor.
## If some energy become available in such circuit, energy interchange starts between
## inductor's magnetic field and capacitor's electric field.
## The interchange process is harmonical oscillation with the oscillation period T.
## Tomson's formula defines this period.
## T = 2 * pi * sqrt(L * C), where
## T is the oscillation period
## L is the inductance of inductor
## C is the capacitance of capacitor
## Pi is 3.1415926538

oscillation_period = Symbol("oscillation_period", units.time)
inductance = Symbol("inductance", units.inductance)
capacitance = Symbol("capacitance", units.capacitance)

law = Eq(oscillation_period, 2 * pi * sqrt(inductance * capacitance))


def print_law() -> str:
    return print_expression(law)


## Derive the same law from the capacitor, charge and self-induction voltage laws

## Let's assume we initially have capacitor charged to U0 voltage. In the zero time we connect this capacitor to inductor in a closed loop.
## So voltage on capacitor is always equals to voltage on inductor and the current through capacitor equals to current through inductor.

#NOTE: this proof is valid for capacitor and inductor in a closed loop without additional voltage source.
#      There are 2 more options: serial connection with external voltage source and parallel connection with external voltage source.
#      Additional proof can be added to show that oscillation period stays the same.

## 1. Prove that capacitor_current(time) = inductor_current(time)

time = Symbol("time", units.time)
capacitor_current = Function("capacitor_current", units.current)
inductor_current = Function("inductor_current", units.current)

local_index_ = Idx("local_index_", (1, 2))
two_currents_law = kirchhoff_law.law.subs(global_index, local_index_).doit()
# capacitor current is in, inductor current is out
two_currents_applied = two_currents_law.subs({
    kirchhoff_law.current[1]: capacitor_current(time),
    kirchhoff_law.current[2]: -1 * inductor_current(time)
})
capacitor_current_applied = solve(two_currents_applied, capacitor_current(time),
    dict=True)[0][capacitor_current(time)]
capacitor_current_eq = Eq(capacitor_current(time), capacitor_current_applied)

assert capacitor_current_eq.lhs == capacitor_current(time)
assert capacitor_current_eq.rhs == inductor_current(time)

## 2. Prove that capacitor_voltage(time) = inductor_voltage(time)

capacitor_voltage = Function("capacitor_voltage", units.voltage)
inductor_voltage = Function("inductor_voltage", units.voltage)

two_voltages_law = kirchhoff_law_2.law.subs(global_index, local_index_).doit()
# capacitor is voltage source, inductor is voltage consumer
two_voltages_applied = two_voltages_law.subs({
    kirchhoff_law_2.voltage[1]: -1 * inductor_voltage(time),
    kirchhoff_law_2.voltage[2]: capacitor_voltage(time)
})
inductor_voltage_applied = solve(two_voltages_applied, inductor_voltage(time),
    dict=True)[0][inductor_voltage(time)]

assert inductor_voltage_applied == capacitor_voltage(time)

## 3. Prove that capacitor current derivative equals to capacitance * (second order derivative of voltage of capacitor)

## charge of capacitor is voltage of capacitor * capacitance
capacitor_charge_law = capacitance_definition.definition.subs({
    capacitance_definition.capacitance: capacitance,
    capacitance_definition.charge: charge_definition.charge(time),
    capacitance_definition.voltage: capacitor_voltage(time)
})
capacitor_charge_applied = solve(capacitor_charge_law, charge_definition.charge(time),
    dict=True)[0][charge_definition.charge(time)]

## I_c(t) = C * U_c'(t)
capacitor_current_law = charge_definition.definition.subs(charge_definition.time, time)
capacitor_current_law = capacitor_current_law.subs(charge_definition.charge(time),
    capacitor_charge_applied)
## I_c'(t) = C * U_c"(t)
capacitor_current_law_derivative = Eq(Derivative(capacitor_current_law.lhs, time),
    Derivative(capacitor_current_law.rhs, time))

## 4. Prove that inductor voltage equals to -1 * capacitance * inductance * (second order derivative of voltage of capacitor)

## Inductor voltage is the self-inductance.
inductor_voltage_law = induction_voltage_definition.definition.subs(
    induction_voltage_definition.time, time)
inductor_voltage_law = inductor_voltage_law.subs({
    induction_voltage_definition.inductance: inductance,
    induction_voltage_definition.self_induction_voltage(time): capacitor_voltage(time),
    induction_voltage_definition.current(time): charge_definition.current(time)
})

derived_law = [inductor_voltage_law, capacitor_current_law_derivative]

## U"(t) = - 1/LC * U(t)
capacitor_voltage_solved = solve(derived_law,
    (Derivative(charge_definition.current(time)), capacitor_voltage(time)),
    dict=True)[0][capacitor_voltage(time)]
voltage_diff_eq = Eq(capacitor_voltage(time), capacitor_voltage_solved)

## 5. Solve differential equation and find period of the harmonic oscillator

## Expected solution for U"(t) = - 1/LC * U(t) is:
## A * e^(i * w * t) + B * e^(-i * w * t), where w = 1 / sqrt(LC)

oscillator_eq = oscillator.definition.subs(oscillator.time, time)
oscillator_eq = oscillator_eq.subs(oscillator.displacement(time), capacitor_voltage(time))
angular_frequency_solved = simplify(
    solve([oscillator_eq, voltage_diff_eq], (oscillator.angular_frequency, capacitor_voltage(time)),
    dict=True)[0][oscillator.angular_frequency])

# 6. Derive period from frequency
period_law = period_definition.law.subs(period_definition.angular_frequency,
    angular_frequency_solved)
period_solved = solve(period_law, period_definition.period, dict=True)[0][period_definition.period]
## Square roots fail to compare with each other. Raise both parts to power of 2 before checking for equality.
assert expr_equals(period_solved**2, law.rhs**2)


@validate_input(inductance_=inductance, capacitance_=capacitance)
@validate_output(oscillation_period)
def calculate_oscillation_period(inductance_: Quantity, capacitance_: Quantity) -> Quantity:
    result_period_expr = solve(law, oscillation_period, dict=True)[0][oscillation_period]
    result_expr = result_period_expr.subs({inductance: inductance_, capacitance: capacitance_})
    return Quantity(result_expr)
