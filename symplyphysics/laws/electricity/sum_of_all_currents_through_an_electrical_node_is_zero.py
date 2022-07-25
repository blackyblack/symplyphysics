from collections import namedtuple
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.laws.electricity import electric_charge_is_constant_in_isolated_system as charge_conservation_law

# Description
## sum(I) = 0
## Where I is a current flowing through electrical node.
## In other words, as electrical charge is being neither created nor accumulated in electrical node,
## sum of all input currents are equal to sum of output current. If we asset input current is positive and output is negative, we gain summary current as 0.
## This property of electrical node is called Kirchhoff law #1.
## Assert there are minimum 2 currents flowing through any node

I_unknown = units.Quantity('I_unknown')
I_known_0 = units.Quantity('I_known_0')
Currents = tuple((I_unknown, I_known_0, ))
law = Eq(sum(Currents),0)

def print():
    return pretty(law, use_unicode=False)

@validate_input(I_known_0 = units.current)
@validate_output(units.current)
def calculate_current(unknown_current_: Quantity) -> Quantity:
    solved = solve(law, Currents.I_unknown, dict=True)[0][Currents.I_unknown]    
    result_expr = solved.subs({Currents.I_unknown : unknown_current_})
    return expr_to_quantity(result_expr, 'unknown_current')
