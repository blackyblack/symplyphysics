"""
Reflection coefficient from load and surge impedance
====================================================

Knowing the load impedance and the surge impedance of the transmission line, it is
possible to calculate the reflection coefficient.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Reflection_coefficient#Relation_to_load_impedance>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

reflection_coefficient = symbols.reflection_coefficient
"""
:symbols:`reflection_coefficient`.
"""

load_impedance = clone_as_symbol(symbols.electrical_impedance, display_symbol="Z_L", display_latex="Z_\\text{L}")
"""
Load :symbols:`electrical_impedance`.
"""

surge_impedance = symbols.surge_impedance
"""
:symbols:`surge_impedance` of the transmission line.
"""

law = Eq(reflection_coefficient,
    (load_impedance - surge_impedance) / (load_impedance + surge_impedance))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(load_impedance_=load_impedance, characteristic_impedance_=surge_impedance)
@validate_output(reflection_coefficient)
def calculate_reflection_coefficient(load_impedance_: Quantity,
    characteristic_impedance_: Quantity) -> float:
    result_expr = solve(law, reflection_coefficient, dict=True)[0][reflection_coefficient]
    result_expr = result_expr.subs({
        load_impedance: load_impedance_,
        surge_impedance: characteristic_impedance_,
    })
    return convert_to_float(result_expr)
