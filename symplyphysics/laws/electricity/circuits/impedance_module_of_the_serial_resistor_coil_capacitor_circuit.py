from sympy import (Eq, solve, sqrt, Idx, expand)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, global_index)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.electricity import capacitor_impedance_from_capacitive_reactance as capacitor_impedance_law
from symplyphysics.laws.electricity import coil_impedance_from_inductive_reactance as coil_impedance_law
from symplyphysics.laws.electricity.circuits import serial_impedance as serial_law

# Description
## Consider an electrical circuit consisting of a capacitor, coil, and resistor connected in series.
## Then you can find the impedance module of such a circuit.

## Law is: |Z| = sqrt(R0^2 + (Xl - Xc)^2), where
## Z - impedance of circuit,
## |Z| - absolute value of Z,
## R0 - resistance of resistor
## Xl - inductive reactance,
## Xc - capacitive reactance.

circuit_resistance = Symbol("circuit_resistance", units.impedance, real=True)

resistance_resistor = Symbol("resistance_resistor", units.impedance, real=True)
capacitive_reactance = Symbol("capacitive_reactance", units.impedance, real=True)
inductive_reactance = Symbol("inductive_reactance", units.impedance, real=True)

law = Eq(circuit_resistance, sqrt(resistance_resistor**2 + (inductive_reactance - capacitive_reactance)**2))

# This law might be derived via "serial_impedance" law, "capacitor_impedance_from_capacitive_reactance" law,
# "coil_impedance_from_inductive_reactance" law.

impedance_law_applied_1 = capacitor_impedance_law.law.subs({
capacitor_impedance_law.capacitive_reactance: capacitive_reactance,
})
capacitive_impedance_derived = solve(impedance_law_applied_1, capacitor_impedance_law.capacitor_impedance, dict=True)[0][capacitor_impedance_law.capacitor_impedance]

impedance_law_applied_2 = coil_impedance_law.law.subs({
coil_impedance_law.inductive_reactance: inductive_reactance,
})
coil_impedance_derived = solve(impedance_law_applied_2, coil_impedance_law.coil_impedance, dict=True)[0][coil_impedance_law.coil_impedance]

local_index = Idx("index_local", (1, 3))
serial_law_applied = serial_law.law.subs(global_index, local_index)
serial_law_applied = serial_law_applied.doit()
circuit_impedance_derived = solve(serial_law_applied, serial_law.serial_impedance, dict=True)[0][serial_law.serial_impedance]
for i, v in enumerate((resistance_resistor, coil_impedance_derived, capacitive_impedance_derived)):
    circuit_impedance_derived = circuit_impedance_derived.subs(serial_law.impedance[i + 1], v)


assert expr_equals(abs(circuit_impedance_derived), expand(law.rhs))


def print_law() -> str:
    return print_expression(law)


@validate_input(resistance_resistor_=resistance_resistor, capacitive_reactance_=capacitive_reactance, inductive_reactance_=inductive_reactance)
@validate_output(circuit_resistance)
def calculate_circuit_impedance_module(resistance_resistor_: Quantity, capacitive_reactance_: Quantity,
    inductive_reactance_: Quantity) -> Quantity:
    result_expr = solve(law, circuit_resistance, dict=True)[0][circuit_resistance]
    result_expr = result_expr.subs({
        resistance_resistor: resistance_resistor_,
        capacitive_reactance: capacitive_reactance_,
        inductive_reactance: inductive_reactance_
    })
    return Quantity(result_expr)
