r"""
Isochoric and isobaric heat capacities of ideal gas
===================================================

Mayer's relation is the relation between heat capacity at constant pressure and that at
constant volume in the case of an ideal gas.

**Notation:**

#. :math:`R` is the molar gas constant.

**Conditions:**

#. Gas is ideal.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

isobaric_heat_capacity = Symbol("isobaric_heat_capacity", units.energy / units.temperature)
"""
Heat capacity of gas at constant pressure.

Symbol:
    :code:`C_p`

Latex:
    :math:`C_p`
"""

isochoric_heat_capacity = Symbol("isochoric_heat_capacity", units.energy / units.temperature)
"""
Heat capacity of gas at constant volume.

Symbol:
    :code:`C_V`

Latex:
    :math:`C_V`
"""

amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)
"""
Amount of gas substance.

Symbol:
    :code:`n`
"""

law = Eq(
    isobaric_heat_capacity - isochoric_heat_capacity,
    amount_of_substance * units.molar_gas_constant,
)
r"""
:code:`C_p - C_V = n * R`

Latex:
    .. math::
        C_p - C_V = n R
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
