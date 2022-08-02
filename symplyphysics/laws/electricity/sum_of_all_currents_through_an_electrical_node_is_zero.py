from collections import namedtuple
from re import I, X
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.laws.electricity import electric_charge_is_constant_in_isolated_system as charge_conservation_law
from sympy import Sum, Function

# Description
## sum(I) = 0
## Where I is a current flowing through electrical node.
## In other words, as electrical charge is being neither created nor accumulated in electrical node,
## sum of all input currents are equal to sum of output current. If we asset input current is positive and output is negative, we gain summary current as 0.
## This property of electrical node is called Kirchhoff law #1.
## Assert there are minimum 2 currents flowing through any node

I0, I1, I2, I3, I4, I5 = symbols('I0, I1, I2, I3, I4, I5')
i = symbols('i', integer = True)
I = Function('I')
def I(x):
    if x == '0':
        return I0
    if x == '1':
        return I1
    if x == '2':
        return I2
    if x == '3':
        return I3
    if x == '4':
        return I4
    if x == '5':
        return I5                    

law = Eq(0, Sum(I(i), (i, 0, 5)).doit())

def print():
    return pretty(law, use_unicode=False)

# @validate_input(I_known_0 = units.current)
# @validate_output(units.current)
# def calculate_current(unknown_current_: Quantity) -> Quantity:
