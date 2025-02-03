from sympy import Eq, solve, exp, acosh
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

attenuation_coefficient = SymbolNew("N", dimensionless)
"""
Attenutation coefficient of the attenuator.

..
    NOTE: why is it dimensionless? maybe it's a different quantity?
"""

first_resistance = clone_as_symbol(symbols.electrical_resistance, subscript="1")
"""
:symbols:`electrical_resistance` of the first section.
"""

second_resistance = clone_as_symbol(symbols.electrical_resistance, subscript="2")
"""
:symbols:`electrical_resistance` of the second section.
"""

law = Eq(attenuation_coefficient, exp(acosh(1 + first_resistance / second_resistance)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(first_resistance_=first_resistance, second_resistance_=second_resistance)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(first_resistance_: Quantity,
    second_resistance_: Quantity) -> float:
    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        first_resistance: first_resistance_,
        second_resistance: second_resistance_,
    })
    return convert_to_float(result_expr)
