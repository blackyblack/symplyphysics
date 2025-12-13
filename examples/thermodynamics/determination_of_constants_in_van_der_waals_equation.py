#!/usr/bin/env python3
r"""
A vessel with a volume of :math:`V = 250 \, \text{mL}` contains :math:`1 \, \text{mol}` of gas.
At temperature :math:`T_1 = 300 \, \text{K}`, the pressure of the gas is
:math:`p_1 = 90 \, \text{atm}`. At temperature :math:`T_2 = 350 \, \text{K}`, the pressure becomes
:math:`p_2 = 110 \, \text{atm}`. Determine the van der Waals constants of this gas.

**Source:** `StudFiles — Example 2.24 <https://studfile.net/preview/1772224/page:8/>`__.
"""

from sympy import solve, Symbol, Eq
from symplyphysics import print_expression, Quantity, units, convert_to
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import (
    van_der_vaals_equation as van_der_waals_law,)
from symplyphysics.reorganized.quantity_relations import (
    quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law,)

temperature_before = Symbol("temperature_before")
temperature_after = Symbol("temperature_after")
volume = Symbol("volume")
amount_of_substance = Symbol("amount_of_substance")
pressure_before = Symbol("pressure_before")
pressure_after = Symbol("pressure_after")

parameter_a = Symbol("parameter_a")
parameter_b = Symbol("parameter_b")

molar_volume = solve(molar_qty_law.law, molar_qty_law.molar_quantity)[0].subs({
    molar_qty_law.extensive_quantity: volume,
    molar_qty_law.amount_of_substance: amount_of_substance,
})

base_variables_of_state_equation = {
    van_der_waals_law.molar_volume: molar_volume,
    van_der_waals_law.attractive_forces_parameter: parameter_a,
    van_der_waals_law.excluded_volume_parameter: parameter_b,
}
base_values_of_state_equation = {
    temperature_before: Quantity(300 * units.kelvins),
    temperature_after: Quantity(350 * units.kelvins),
    volume: Quantity(0.25 * units.liters),
    amount_of_substance: Quantity(1 * units.mole),
    pressure_before: Quantity(90 * units.atmospheres),
    pressure_after: Quantity(110 * units.atmospheres),
}

state_equation_before = van_der_waals_law.law.subs({
    van_der_waals_law.temperature: temperature_before,
    van_der_waals_law.pressure: pressure_before,
}).subs(base_variables_of_state_equation)

state_equation_after = van_der_waals_law.law.subs({
    van_der_waals_law.temperature: temperature_after,
    van_der_waals_law.pressure: pressure_after,
}).subs(base_variables_of_state_equation)

solved = solve(
    (state_equation_before, state_equation_after),
    (parameter_a, parameter_b),
    dict=True,
)[0]
parameter_a_value = solved[parameter_a]
parameter_b_value = solved[parameter_b]

answer1 = Eq(parameter_a, parameter_a_value)
answer2 = Eq(parameter_b, parameter_b_value)
print(f"Equation for parameter a:\n{print_expression(answer1)}")
print(f"Equation for parameter b:\n{print_expression(answer2)}")

parameter_a_ = answer1.subs(base_values_of_state_equation).rhs
parameter_b_ = answer2.subs(base_values_of_state_equation).rhs

answer1_value = convert_to(
    Quantity(parameter_a_),
    units.pascals * (units.meters**3 / units.mole)**2,
)
answer2_value = convert_to(
    Quantity(parameter_b_),
    units.liters / units.mole,
)
print(f"Parameter a is: {answer1_value} Pa * (m^3 / mol)^2")
print(f"Parameter b is: {answer2_value} L / mol")
