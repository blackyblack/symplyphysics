from sympy import (Eq, solve, sqrt,)
from symplyphysics import (units, Quantity, Symbol, Vector, print_expression, validate_input,
    validate_output, vector_magnitude,)
from symplyphysics.core.expr_comparisons import expr_equals

# Description
## Consider an electrical circuit consisting of a capacitor, coil, and resistor connected in series.
## Then you can find the total resistance of such a circuit.

## Law is: R = sqrt(R0^2 + (Xl - Xc)^2), where
## R - resistance of circuit,
## R0 - resistance of resistor
## Xl - inductive resistance,
## Xc - capacitive resistance.

circuit_resistance = Symbol("circuit_resistance", units.impedance)

resistance_resistor = Symbol("resistance_resistor", units.impedance)
capacitive_resistance = Symbol("capacitive_resistance", units.impedance)
inductive_resistance = Symbol("inductive_resistance", units.impedance)

law = Eq(circuit_resistance, sqrt(resistance_resistor**2 + (inductive_resistance - capacitive_resistance)**2))

# This law might be derived via impedance definition. The impedance is Z = R0 + j * X.
# Then the resistance of the circuit is equal to the impedance modulus and is equal to sqrt(R0^2 + X^2).
# X = Xc - Xl
# Then, as a result, the resistance of the circuit is equal to sqrt(R0^2 + (Xl - Xc)^2).

# The complex impedance can be represented as a vector with two components: the real and imaginary parts.
# Let's find the modulus of this vector and get the resistance of the circuit.

resistance_vector = Vector([resistance_resistor, inductive_resistance - capacitive_resistance])
resistance_derived = vector_magnitude(resistance_vector)

# Check if derived resistance is same as declared.
assert expr_equals(resistance_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(resistance_resistor_=resistance_resistor, capacitive_resistance_=capacitive_resistance, inductive_resistance_=inductive_resistance)
@validate_output(circuit_resistance)
def calculate_circuit_resistance(resistance_resistor_: Quantity, capacitive_resistance_: Quantity,
    inductive_resistance_: Quantity) -> Quantity:
    result_expr = solve(law, circuit_resistance, dict=True)[0][circuit_resistance]
    result_expr = result_expr.subs({
        resistance_resistor: resistance_resistor_,
        capacitive_resistance: capacitive_resistance_,
        inductive_resistance: inductive_resistance_
    })
    return Quantity(result_expr)
