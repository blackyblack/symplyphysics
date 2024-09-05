from sympy import (Derivative, Eq, Idx, solve, exp, simplify)
from symplyphysics import (units, Quantity, Symbol, Function, print_expression, validate_input,
    validate_output, global_index)
from symplyphysics.definitions import current_is_charge_derivative as charge_definition
from symplyphysics.definitions import capacitance_from_charge_and_voltage as capacitance_definition
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohms_law
from symplyphysics.laws.electricity.circuits import sum_of_currents_through_junction_is_zero as kirchhoff_law
from symplyphysics.laws.electricity.circuits import sum_of_voltages_in_loop_is_zero as kirchhoff_law_2

# Description
## RC integrator is a circuit with capacitor and resistor in series. Initial_voltage is applied to whole circuit and integrated voltage is obtained from capacitor.
## If some voltage is applied to RC integrator, capacitor starts to charge and it's voltage rises.
## Capacitor is charging, but it's voltage will never reach initial voltage.

# Law: Uc(t) = U0 * (1 - e^(-t / (R * C)))
## Where:
## t - time
## e - exponent
## Uc - capacitor voltage
## R - resistor impedance
## C - capacitor capacitance

time = Symbol("time", units.time)
initial_voltage = Symbol("initial_voltage", units.voltage)
capacitance = Symbol("capacitance", units.capacitance)
resistance = Symbol("resistance", units.impedance)
capacitor_voltage = Function("capacitor_voltage", units.voltage)

law = Eq(capacitor_voltage(time), initial_voltage * (1 - exp(-time / (capacitance * resistance))))

## Derive the same law from the Ohms and Kirchhoff laws

capacitor_current = Function("capacitor_current", units.current)
resistor_current = Function("resistor_current", units.current)
resistor_voltage = Function("resistor_voltage", units.voltage)

local_index_ = Idx("local_index_", (1, 2))
two_currents_law = kirchhoff_law.law.subs(global_index, local_index_).doit()
# capacitor current is in, resistor current is out
two_currents_applied = two_currents_law.subs({
    kirchhoff_law.current[1]: capacitor_current(time),
    kirchhoff_law.current[2]: -1 * resistor_current(time)
})
capacitor_current_applied = solve(two_currents_applied, capacitor_current(time),
    dict=True)[0][capacitor_current(time)]
capacitor_current_eq = Eq(capacitor_current(time), capacitor_current_applied)

## 1. Prove that capacitor_current(time) = resistor_current(time)
assert capacitor_current_eq.lhs == capacitor_current(time)
assert capacitor_current_eq.rhs == resistor_current(time)

local_index_ = Idx("local_index_", (1, 3))
three_voltages_law = kirchhoff_law_2.law.subs(global_index, local_index_).doit()
# initial_voltage is voltage source, capacitor and resistor are voltage consumers
three_voltages_applied = three_voltages_law.subs({
    kirchhoff_law_2.voltage[1]: -1 * capacitor_voltage(time),
    kirchhoff_law_2.voltage[2]: -1 * resistor_voltage(time),
    kirchhoff_law_2.voltage[3]: initial_voltage
})
resistor_voltage_applied = solve(three_voltages_applied, resistor_voltage(time),
    dict=True)[0][resistor_voltage(time)]
resistor_voltage_eq = Eq(resistor_voltage(time), resistor_voltage_applied)

## 2. Prove that resistor_voltage(time) = initial_voltage - capacitor_voltage(time)
assert resistor_voltage_eq.rhs == initial_voltage - capacitor_voltage(time)

# use resistor_voltage as proven in resistor_voltage_eq
# use charge_definition.current since it is same on resistor and capacitor as proven in capacitor_current_eq
resistor_ohm_eq = ohms_law.law.subs({
    ohms_law.voltage: initial_voltage - capacitor_voltage(time),
    ohms_law.resistance: resistance,
    ohms_law.current: charge_definition.current(time)
})
capacitance_eq = capacitance_definition.definition.subs({
    capacitance_definition.capacitance: capacitance,
    capacitance_definition.charge: charge_definition.charge(time),
    capacitance_definition.voltage: capacitor_voltage(time)
})
charge_eq = charge_definition.definition.subs(charge_definition.time, time)

derived_law = [resistor_ohm_eq, capacitance_eq, charge_eq]

solved_charge_function = solve(derived_law,
    (capacitor_voltage(time), charge_definition.current(time), charge_definition.charge(time)),
    dict=True)[0][charge_definition.charge(time)]
charge_diff_eq = Eq(charge_definition.charge(time), solved_charge_function)

## 3. Prove that charge(time) = capacitance * initial_voltage - capacitance * resistance * Derivative(charge(time), time))
## Q(t) = U0 * C - R * C * dQ(t) / dt
capacitor_charge_function = initial_voltage * capacitance - resistance * capacitance * Derivative(
    charge_definition.charge(time), time)
assert simplify(charge_diff_eq.rhs - capacitor_charge_function) == 0

## 4. Convert charge to capacitor voltage
capacitor_voltage_solved = solve(capacitance_eq, charge_definition.charge(time),
    dict=True)[0][charge_definition.charge(time)]
voltage_diff_eq = charge_diff_eq.subs(charge_definition.charge(time), capacitor_voltage_solved)

## 5. Solve differential equation
# HACK: use known solution since sympy.dsolve() gives us another result
voltage_diff_solution = voltage_diff_eq.subs(capacitor_voltage(time), law.rhs)
assert simplify(voltage_diff_solution.lhs - voltage_diff_solution.rhs) == 0


def print_law() -> str:
    return print_expression(law)


@validate_input(initial_voltage_=initial_voltage,
    capacitance_=capacitance,
    resistance_=resistance,
    time_=time)
@validate_output(capacitor_voltage)
def calculate_capacitor_voltage(initial_voltage_: Quantity, capacitance_: Quantity,
    resistance_: Quantity, time_: Quantity) -> Quantity:
    capacitor_voltage_expr = solve(law, capacitor_voltage(time),
        dict=True)[0][capacitor_voltage(time)]
    result_expr = capacitor_voltage_expr.subs({
        initial_voltage: initial_voltage_,
        resistance: resistance_,
        capacitance: capacitance_,
        time: time_
    })
    return Quantity(result_expr)
