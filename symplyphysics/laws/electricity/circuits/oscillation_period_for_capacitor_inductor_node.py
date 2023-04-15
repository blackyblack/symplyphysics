from sympy import sqrt, pi, cos, Derivative
from symplyphysics import (
    symbols, Eq, pretty, solve, simplify, Quantity, units, Function,
    validate_input, validate_output, expr_to_quantity
)

from symplyphysics.laws.electricity.circuits import sum_of_all_currents_through_an_electrical_node_is_zero as kirchhoff_law
from symplyphysics.laws.electricity.circuits import sum_of_all_voltages_in_loop_is_zero as kirchhoff_law_2
from symplyphysics.definitions import capacitance_from_charge_and_voltage as capacitance_definition
from symplyphysics.definitions import self_induction_voltage_is_current_derivative as induction_voltage_definition
from symplyphysics.definitions import current_is_charge_derivative as charge_definition
from symplyphysics.definitions import period_from_circular_frequency as period_definition

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

oscillation_period, inductance, capacitance = symbols('oscillation_period inductance capacitance')
law = Eq(oscillation_period, 2 * pi * sqrt(inductance * capacitance))

def print():
    return pretty(law, use_unicode=False)

## Derive the same law from the capacitor, charge and self-induction voltage laws

## Let's assume we initially have capacitor charged to U0 voltage. In the zero time we connect this capacitor to inductor in a closed loop.
## So voltage on capacitor is always equals to voltage on inductor and the current through capacitor equals to current through inductor.

#NOTE: this proof is valid for capacitor and inductor in a closed loop without additional voltage source.
#      There are 2 more options: serial connection with external voltage source and parallel connection with external voltage source.
#      Additional proof can be added to show that oscillation period stays the same.

## 1. Prove that capacitor_current_function(time) = inductor_current_function(time)

time = symbols('time')
capacitor_current_function, inductor_current_function = symbols('capacitor_current_function inductor_current_function', cls = Function)

two_currents_law = kirchhoff_law.law.subs(kirchhoff_law.currents_total, 2).doit()
# capacitor current is in, inductor current is out
two_currents_applied = two_currents_law.subs({
    kirchhoff_law.current[1]: capacitor_current_function(time),
    kirchhoff_law.current[2]: -1 * inductor_current_function(time)})
capacitor_current_applied = solve(two_currents_applied, capacitor_current_function(time), dict=True)[0][capacitor_current_function(time)]
capacitor_current_eq = Eq(capacitor_current_function(time), capacitor_current_applied)

assert capacitor_current_eq.lhs == capacitor_current_function(time)
assert capacitor_current_eq.rhs == inductor_current_function(time)

## 2. Prove that capacitor_voltage_function(time) = inductor_voltage_function(time)

capacitor_voltage_function = symbols('capacitor_voltage_function', cls = Function)
inductor_voltage_function = symbols('inductor_voltage_function', cls = Function)

two_voltages_law = kirchhoff_law_2.law.subs(kirchhoff_law_2.voltages_total, 2).doit()
# capacitor is voltage source, inductor is voltage consumer
two_voltages_applied = two_voltages_law.subs({
    kirchhoff_law_2.voltage[1]: -1 * inductor_voltage_function(time),
    kirchhoff_law_2.voltage[2]: capacitor_voltage_function(time)})
inductor_voltage_applied = solve(two_voltages_applied, inductor_voltage_function(time), dict=True)[0][inductor_voltage_function(time)]

assert inductor_voltage_applied == capacitor_voltage_function(time)

## 3. Prove that capacitor current derivative equals to capacitance * (second order derivative of voltage of capacitor)

## charge of capacitor is voltage of capacitor * capacitance
capacitor_charge_law = capacitance_definition.definition.subs(
    {capacitance_definition.capacitance: capacitance,
    capacitance_definition.charge: charge_definition.charge(time),
    capacitance_definition.voltage: capacitor_voltage_function(time)})
capacitor_charge_applied = solve(capacitor_charge_law, charge_definition.charge(time), dict=True)[0][charge_definition.charge(time)]

## I_c(t) = C * U_c'(t)
capacitor_current_law = charge_definition.definition.subs(charge_definition.time, time)
capacitor_current_law = capacitor_current_law.subs(charge_definition.charge(time), capacitor_charge_applied)
## I_c'(t) = C * U_c"(t)
capacitor_current_law_derivative = Eq(Derivative(capacitor_current_law.lhs, time), Derivative(capacitor_current_law.rhs, time))

## 4. Prove that inductor voltage equals to -1 * capacitance * inductance * (second order derivative of voltage of capacitor)

## Inductor voltage is the self-inductance.
inductor_voltage_law = induction_voltage_definition.definition.subs(induction_voltage_definition.time, time)
inductor_voltage_law = inductor_voltage_law.subs(
    {induction_voltage_definition.inductance: inductance,
    induction_voltage_definition.self_induction_voltage(time): capacitor_voltage_function(time),
    induction_voltage_definition.current(time): charge_definition.current(time)})

derived_law = [inductor_voltage_law, capacitor_current_law_derivative]

## U"(t) = - 1/LC * U(t)
capacitor_voltage_solved = solve(derived_law, (Derivative(charge_definition.current(time)), capacitor_voltage_function(time)), dict=True)[0][capacitor_voltage_function(time)]
voltage_diff_eq = Eq(capacitor_voltage_function(time), capacitor_voltage_solved)

## 5. Solve differential equation and find period of the harmonic oscillator

## Expected solution for U"(t) = - 1/LC * U(t) is:
## A * e^(i * w * t) + B * e^(-i * w * t), where w = 1 / sqrt(LC)
## This form is known to be represented as cosine function:
## A * cos(w * t + phi), or sine function:
## A * sin(w * t + phi + pi / 2)

## Let's check this solution:
amplitude, phase = symbols('amplitude phase')
frequency = 1 / sqrt(inductance * capacitance)
voltage_diff_solution = amplitude * cos(frequency * time + phase)
voltage_diff_solved = voltage_diff_eq.subs({capacitor_voltage_function(time): voltage_diff_solution})

assert simplify(voltage_diff_solved.lhs - voltage_diff_solved.rhs) == 0

# 6. Derive period from frequency
period_law = period_definition.definition.subs(period_definition.circular_frequency, frequency)
period_solved = solve(period_law, period_definition.period, dict=True)[0][period_definition.period]

assert simplify(period_solved - law.rhs) == 0


@validate_input(inductance_=units.inductance, capacitance_=units.capacitance)
@validate_output(units.time)
def calculate_oscillation_period(inductance_: Quantity, capacitance_: Quantity) -> Quantity:
    result_period_expr = solve(law, oscillation_period, dict=True)[0][oscillation_period]
    result_expr = result_period_expr.subs({inductance: inductance_, capacitance: capacitance_})
    return expr_to_quantity(result_expr, 'oscillation_period')
