"""
Attenuation of three link microwave attenuator
==============================================

Microwave attenuators are used to attenuate the microwave signal. For a three-link
T-type attenuator or a π-type attenuator, the signal attenuation coefficient is
calculated from the ratio of resistances of the resistors.

..
    TODO: find link
    TODO: fix file name
"""

from sympy import Eq, solve, exp, acosh
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

attenuation = symbols.attenuation
"""
:symbols:`attenuation` of the attenuator.
"""

first_resistance = clone_as_symbol(symbols.electrical_resistance, subscript="1")
"""
:symbols:`electrical_resistance` of the first section.
"""

second_resistance = clone_as_symbol(symbols.electrical_resistance, subscript="2")
"""
:symbols:`electrical_resistance` of the second section.
"""

law = Eq(attenuation, exp(acosh(1 + first_resistance / second_resistance)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(first_resistance_=first_resistance, second_resistance_=second_resistance)
@validate_output(attenuation)
def calculate_attenuation_coefficient(first_resistance_: Quantity,
    second_resistance_: Quantity) -> float:
    result_expr = solve(law, attenuation, dict=True)[0][attenuation]
    result_expr = result_expr.subs({
        first_resistance: first_resistance_,
        second_resistance: second_resistance_,
    })
    return convert_to_float(result_expr)
