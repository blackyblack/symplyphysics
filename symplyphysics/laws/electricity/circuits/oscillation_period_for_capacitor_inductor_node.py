from sympy import sqrt, pi
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

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

## Let's assume we initially have capacitor charged to U0 voltage. In the zero time we connect this capacitor to inductor.
## So voltage on capacitor is always equals to voltage on inductor and the current through capacitor equals to current through inductor.

## Current through capacitor is charge derivative, charge of capacitor is voltage of capacitor * capacitance
## I_c(t) = C * U_c'(t)
## Derive rsh and lhs:
## I_c'(t) = C * U_c"(t)

## Inductor voltage is the self-inductance.
## U_i(t) = - L * I_i'(t) => I_i'(t) = U_i(t) / L

## U-C = U_l = U, I_c = I_l = I

## C * U"(t) = -L U(t). U"(t) = - 1/LC * U(t)
## U(t) = U0 * cos(wt), where w = 1/sqrt(LC) and the period T is 2*pi/w = 2 * pi * sqrt(LC)

@validate_input(inductance_=units.inductance, capacitance_=units.capacitance)
@validate_output(units.time)
def calculate_oscillation_period(inductance_: Quantity, capacitance_: Quantity) -> Quantity:
    result_period_expr = solve(law, oscillation_period, dict=True)[0][oscillation_period]
    result_expr = result_period_expr.subs({inductance: inductance_, capacitance: capacitance_})
    return expr_to_quantity(result_expr, 'oscillation_period')
