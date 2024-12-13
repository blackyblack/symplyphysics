from sympy import (Eq, solve, sqrt, Idx, expand)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, global_index)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.electricity.circuits import capacitor_impedance_from_capacitive_reactance as capacitor_impedance_law
from symplyphysics.laws.electricity.circuits import coil_impedance_from_inductive_reactance as coil_impedance_law
from symplyphysics.laws.electricity.circuits import impedance_in_serial_connection as serial_law

# Description
## Consider an electrical circuit consisting of a capacitor, coil, and resistor connected in series.
## Then you can find the impedance module of such a circuit.

## Law is: |Z| = sqrt(R0^2 + (Xl - Xc)^2), where
## Z - impedance of circuit,
## |Z| - absolute value of Z,
## R0 - resistance of resistor
## Xl - inductive reactance,
## Xc - capacitive reactance.

# Links: derivable from first formula <https://en.wikipedia.org/wiki/Electrical_impedance#Resistance_vs_reactance>

circuit_resistance = Symbol("circuit_resistance", units.impedance, real=True)

resistance_resistor = Symbol("resistance_resistor", units.impedance, real=True)
capacitive_reactance = Symbol("capacitive_reactance", units.impedance, real=True)
inductive_reactance = Symbol("inductive_reactance", units.impedance, real=True)

law = Eq(circuit_resistance,
    sqrt(resistance_resistor**2 + (inductive_reactance - capacitive_reactance)**2))

# This law might be derived via "serial_impedance" law, "capacitor_impedance_from_capacitive_reactance" law,
# "coil_impedance_from_inductive_reactance" law.

_impedance_law_applied_1 = capacitor_impedance_law.law.subs({
    capacitor_impedance_law.capacitive_reactance: capacitive_reactance,
})
_capacitive_impedance_derived = solve(_impedance_law_applied_1,
    capacitor_impedance_law.capacitor_impedance,
    dict=True)[0][capacitor_impedance_law.capacitor_impedance]

_impedance_law_applied_2 = coil_impedance_law.law.subs({
    coil_impedance_law.inductive_reactance: inductive_reactance,
})
_coil_impedance_derived = solve(_impedance_law_applied_2,
    coil_impedance_law.coil_impedance,
    dict=True)[0][coil_impedance_law.coil_impedance]

_local_index = Idx("index_local", (1, 3))
_serial_law_applied = serial_law.law.subs(global_index, _local_index)
_serial_law_applied = _serial_law_applied.doit()
_circuit_impedance_derived = solve(_serial_law_applied, serial_law.total_impedance,
    dict=True)[0][serial_law.total_impedance]
for i, v in enumerate(
    (resistance_resistor, _coil_impedance_derived, _capacitive_impedance_derived)):
    _circuit_impedance_derived = _circuit_impedance_derived.subs(serial_law.impedance[i + 1], v)

assert expr_equals(abs(_circuit_impedance_derived), expand(law.rhs))


@validate_input(resistance_resistor_=resistance_resistor,
    capacitive_reactance_=capacitive_reactance,
    inductive_reactance_=inductive_reactance)
@validate_output(circuit_resistance)
def calculate_circuit_impedance_module(resistance_resistor_: Quantity,
    capacitive_reactance_: Quantity, inductive_reactance_: Quantity) -> Quantity:
    result_expr = solve(law, circuit_resistance, dict=True)[0][circuit_resistance]
    result_expr = result_expr.subs({
        resistance_resistor: resistance_resistor_,
        capacitive_reactance: capacitive_reactance_,
        inductive_reactance: inductive_reactance_
    })
    return Quantity(result_expr)
