"""
Isochoric and isobaric heat capacities of ideal gas
===================================================

The **Mayer's relation** is the relation between heat capacity at constant pressure and that at
constant volume in the case of an ideal gas.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Conditions:**

#. Gas is ideal.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    quantities,
    symbols,
    clone_as_symbol,
)

isobaric_heat_capacity = clone_as_symbol(symbols.heat_capacity, subscript="p")
"""
:symbols:`heat_capacity` of gas at constant :symbols:`pressure`.
"""

isochoric_heat_capacity = clone_as_symbol(symbols.heat_capacity, subscript="V")
"""
:symbols:`heat_capacity` of gas at constant :symbols:`volume`.
"""

amount_of_substance = symbols.amount_of_substance
"""
:symbols:`amount_of_substance`.
"""

law = Eq(
    isobaric_heat_capacity - isochoric_heat_capacity,
    amount_of_substance * quantities.molar_gas_constant,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    isobaric_heat_capacity_=isobaric_heat_capacity,
    amount_of_substance_=amount_of_substance,
)
@validate_output(isochoric_heat_capacity)
def calculate_isochoric_heat_capacity(
    isobaric_heat_capacity_: Quantity,
    amount_of_substance_: Quantity,
) -> Quantity:
    expr = solve(law, isochoric_heat_capacity)[0]
    result = expr.subs({
        isobaric_heat_capacity: isobaric_heat_capacity_,
        amount_of_substance: amount_of_substance_,
    })
    return Quantity(result)
