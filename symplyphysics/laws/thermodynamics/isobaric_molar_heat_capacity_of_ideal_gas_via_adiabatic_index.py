"""
Isobaric molar heat capacity of ideal gas via adiabatic index
=============================================================

For ideal gases the isobaric and isochoric heat capacities can be calculated via the molar gas
constant and the adiabatic index of the gas.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Conditions:**

#. The gas is ideal.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Heat_capacity_ratio#Ideal-gas_relations>`__.
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
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import heat_capacity_ratio
from symplyphysics.laws.thermodynamics import (
    isochoric_and_isobaric_heat_capacity_of_ideal_gas as mayers_law,)
from symplyphysics.laws.quantities import (
    quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law,)

isobaric_molar_heat_capacity = clone_as_symbol(
    symbols.molar_heat_capacity,
    display_symbol="c_pm",
    display_latex="c_{p, \\text{m}}",
)
"""
:symbols:`molar_heat_capacity` at constant :symbols:`pressure`.
"""

adiabatic_index = symbols.adiabatic_index
"""
:symbols:`adiabatic_index`, or :doc:`heat capacity ratio <definitions.heat_capacity_ratio>`
of the gas.
"""

law = Eq(isobaric_molar_heat_capacity,
    quantities.molar_gas_constant * adiabatic_index / (adiabatic_index - 1))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from the Mayer's law for heat capacities and the definition of adiabatic index

_isobaric_heat = mayers_law.isobaric_heat_capacity
_isochoric_heat = mayers_law.isochoric_heat_capacity

_adiabatic_index_eqn = heat_capacity_ratio.definition.subs({
    heat_capacity_ratio.isobaric_heat_capacity: _isobaric_heat,
    heat_capacity_ratio.isochoric_heat_capacity: _isochoric_heat,
    heat_capacity_ratio.heat_capacity_ratio: adiabatic_index,
})

_isobaric_heat_expr = solve(
    (mayers_law.law, _adiabatic_index_eqn),
    (_isobaric_heat, _isochoric_heat),
    dict=True,
)[0][_isobaric_heat]

_molar_isobaric_heat = solve(molar_qty_law.law, molar_qty_law.molar_quantity)[0].subs({
    molar_qty_law.extensive_quantity: _isobaric_heat_expr,
    molar_qty_law.amount_of_substance: mayers_law.amount_of_substance,
})

assert expr_equals(_molar_isobaric_heat, law.rhs)


@validate_input(adiabatic_index_=adiabatic_index)
@validate_output(isobaric_molar_heat_capacity)
def calculate_isobaric_molar_heat_capacity(adiabatic_index_: float) -> Quantity:
    result = law.rhs.subs(adiabatic_index, adiabatic_index_)
    return Quantity(result)
