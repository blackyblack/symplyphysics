from sympy import (Eq, solve, sqrt, I)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output,)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.definitions import impedance_is_resistance_and_reactance as impedance_law
from symplyphysics.laws.electricity.circuits import resistivity_of_serial_resistors as serial_law

# Description
## Consider an electrical circuit consisting of a capacitor, coil, and resistor connected in series.
## Then you can find the total resistance of such a circuit.

## Law is: R = sqrt(R0^2 + (Xl - Xc)^2), where
## R - resistance of circuit,
## R0 - resistance of resistor
## Xl - inductive reactance,
## Xc - capacitive reactance.

circuit_resistance = Symbol("circuit_resistance", units.impedance, real=True)

resistance_resistor = Symbol("resistance_resistor", units.impedance, real=True)
capacitive_reactance = Symbol("capacitive_reactance", units.impedance, real=True)
inductive_reactance = Symbol("inductive_reactance", units.impedance, real=True)

law = Eq(circuit_resistance, sqrt(resistance_resistor**2 + (inductive_reactance - capacitive_reactance)**2))

# This law might be derived via impedance definition and "resistivity_of_serial_resistors" law.

impedance_law_applied_1 = impedance_law.definition.subs({
impedance_law.impedance: -I * capacitive_reactance,
impedance_law.resistance: 0,
})
capacitive_reactance_derived = solve(impedance_law_applied_1, impedance_law.reactance, dict=True)[0][impedance_law.reactance]

impedance_law_applied_2 = impedance_law.definition.subs({
impedance_law.impedance: I * inductive_reactance,
impedance_law.resistance: 0,
})
inductive_reactance_derived = solve(impedance_law_applied_2, impedance_law.reactance, dict=True)[0][impedance_law.reactance]

serial_law_applied = serial_law.law.subs({
    serial_law.resistances: (inductive_reactance_derived, capacitive_reactance_derived)
})
circuit_reactance_derived = solve(serial_law_applied, serial_law.serial_resistance, dict=True)[0][serial_law.serial_resistance]

impedance_law_applied = impedance_law.definition.subs({
impedance_law.resistance: resistance_resistor,
impedance_law.reactance: circuit_reactance_derived,
})
circuit_resistance_derived = solve(impedance_law_applied, impedance_law.impedance, dict=True)[0][impedance_law.impedance]

assert expr_equals(abs(circuit_resistance_derived), law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(resistance_resistor_=resistance_resistor, capacitive_reactance_=capacitive_reactance, inductive_reactance_=inductive_reactance)
@validate_output(circuit_resistance)
def calculate_circuit_resistance(resistance_resistor_: Quantity, capacitive_reactance_: Quantity,
    inductive_reactance_: Quantity) -> Quantity:
    result_expr = solve(law, circuit_resistance, dict=True)[0][circuit_resistance]
    result_expr = result_expr.subs({
        resistance_resistor: resistance_resistor_,
        capacitive_reactance: capacitive_reactance_,
        inductive_reactance: inductive_reactance_
    })
    return Quantity(result_expr)
