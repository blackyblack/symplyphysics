
from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)


# Description
## Ohm's law for the full section of the circuit for direct current is formulated: 
## The current in the circuit is directly proportional to the sum of the EMF of the circuit and inversely
## proportional to the sum of the resistances of the source and the circuit elements
## law: I = EMF / (R + r)
## Where I is the current flowing through circuit 
## EMF - is electromotive force
## R is impedance of this circuit
## r is impedance of the power supply

current = Symbol("current", units.current)
EMF = Symbol("EMF", units.voltage)
resistance = Symbol("resistance", units.impedance)
res_power = Symbol("resistance", units.impedance)

law = Eq(current, EMF / (resistance + res_power))

def print_law() -> str:
    return print_expression(law)


@validate_input(EMF_=EMF, resistance_=resistance, res_power_=res_power)
@validate_output(current)
def calculate_current(EMF_: EMF, resistance_: resistance, res_power_: res_power) -> Quantity:
    result_current_expr = solve(law, current, dict=True)[0][current]
    result_expr = result_current_expr.subs({EMF: EMF_, resistance: resistance_, res_power: res_power_})
    return Quantity(result_expr)
