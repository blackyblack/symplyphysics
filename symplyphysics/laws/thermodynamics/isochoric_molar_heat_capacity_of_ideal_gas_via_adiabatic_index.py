r"""
Isochoric molar heat capacity of ideal gas via adiabatic index
==============================================================

For ideal gases the isobaric and isochoric heat capacities can be calculated via the molar gas
constant and the adiabatic index of the gas.

**Notation:**

#. :math:`R` is the molar gas constant.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import heat_capacity_ratio
from symplyphysics.laws.thermodynamics import (
    isochoric_and_isobaric_heat_capacity_of_ideal_gas as mayers_law,)
from symplyphysics.laws.quantities import (
    quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law,)

isochoric_molar_heat_capacity = Symbol(
    "isochoric_molar_heat_capacity", units.energy / (units.temperature * units.amount_of_substance))
r"""
Heat capacity of ideal gas at constant volume per unit amount of substance.

Symbol:
    :code:`C_V`

Latex:
    :math:`C_V`
"""

adiabatic_index = Symbol("adiabatic_index", dimensionless)
r"""
Adiabatic index, also known as :doc:`heat capacity ratio <definitions.heat_capacity_ratio>`, of the gas.

Symbol:
    :code:`gamma`

Latex:
    :math:`\gamma`
"""

law = Eq(isochoric_molar_heat_capacity, units.molar_gas_constant / (adiabatic_index - 1))
r"""
:code:`C_V = R / (gamma - 1)`

Latex:
    .. math::
        C_V = \frac{R}{\gamma - 1}
"""

# Derive law from the Mayer's law for heat capacities and the definition of adiabatic index

_isobaric_heat = mayers_law.isobaric_heat_capacity
_isochoric_heat = mayers_law.isochoric_heat_capacity

_adiabatic_index_eqn = heat_capacity_ratio.definition.subs({
    heat_capacity_ratio.isobaric_heat_capacity: _isobaric_heat,
    heat_capacity_ratio.isochoric_heat_capacity: _isochoric_heat,
    heat_capacity_ratio.heat_capacity_ratio: adiabatic_index,
})

_isochoric_heat_expr = solve(
    (mayers_law.law, _adiabatic_index_eqn),
    (_isobaric_heat, _isochoric_heat),
    dict=True,
)[0][_isochoric_heat]

_molar_isochoric_heat = solve(molar_qty_law.law, molar_qty_law.molar_quantity)[0].subs({
    molar_qty_law.extensive_quantity: _isochoric_heat_expr,
    molar_qty_law.amount_of_substance: mayers_law.amount_of_substance,
})

assert expr_equals(_molar_isochoric_heat, law.rhs)


@validate_input(adiabatic_index_=adiabatic_index)
@validate_output(isochoric_molar_heat_capacity)
def calculate_isochoric_molar_heat_capacity(adiabatic_index_: float) -> Quantity:
    result = law.rhs.subs(adiabatic_index, adiabatic_index_)
    return Quantity(result)
