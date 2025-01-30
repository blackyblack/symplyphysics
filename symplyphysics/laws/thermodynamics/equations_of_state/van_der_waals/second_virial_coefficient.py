r"""
Second virial coefficient
=========================

The *second virial coefficient* is a coefficient appearing in the virial equation of state of
a gas which describes the deviation from ideal gas behaviour. This coefficient describes pair
interaction between molecules of the substance and can be found if the pair interaction potential
is known, but it can also be derived from the equation of state as a series expansion with respect
to inverse molar volume or equivalently, to molar density.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Notes:**

#. Also see :ref:`Virial equation`.

**Conditions:**

#. The gas density is small enough within the context of perturbation theory.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Van_der_Waals_equation#Virial_expansion>`__.
"""

from sympy import Eq, solve, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
    SymbolNew,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import compressibility_factor_is_deviation_from_ideal_gas as compressibility_def
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import equation as vdw_eqn
from symplyphysics.laws.quantities import (
    quantity_is_volumetric_density_times_volume as density_qty_law,
    quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law,
)

second_virial_coefficient = SymbolNew("C_2", units.volume / units.amount_of_substance)
"""
Second virial coefficient of the virial expansion.
"""

attractive_forces_parameter = symbols.attractive_forces_parameter
"""
:symbols:`attractive_forces_parameter`.
"""

excluded_volume_parameter = symbols.excluded_volume_parameter
"""
:symbols:`excluded_volume_parameter`.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the van der Waals fluid.
"""

law = Eq(
    second_virial_coefficient, excluded_volume_parameter - attractive_forces_parameter /
    (quantities.molar_gas_constant * temperature))
"""
:laws:symbol::

:laws:latex::
"""

# Derive from the van der Waals equation of state and the virial equation

_volume = compressibility_def.volume
_mole_count = compressibility_def.amount_of_substance
_molar_density = SymSymbol("molar_density")

_molar_volume = solve(molar_qty_law.law, molar_qty_law.molar_quantity)[0].subs({
    molar_qty_law.extensive_quantity: _volume,
    molar_qty_law.amount_of_substance: _mole_count,
})

_pressure_expr = solve(vdw_eqn.law, vdw_eqn.pressure)[0].subs({
    vdw_eqn.molar_volume: _molar_volume,
    vdw_eqn.temperature: temperature,
    vdw_eqn.attractive_forces_parameter: attractive_forces_parameter,
    vdw_eqn.excluded_volume_parameter: excluded_volume_parameter,
})

_compressibility_via_volume = compressibility_def.definition.rhs.subs({
    compressibility_def.pressure: _pressure_expr,
    compressibility_def.temperature: temperature,
})

_mole_count_expr = density_qty_law.law.rhs.subs({
    density_qty_law.volumetric_density: _molar_density,
    density_qty_law.volume: _volume,
})

_compressibility_via_density = _compressibility_via_volume.subs(_mole_count,
    _mole_count_expr).simplify()

# Virial equation is a power series of `Z` around `rho = 0`, therefore we can use
# the formula for Taylor series coefficients to get the coefficient at `rho**2`.
_second_virial_coefficient = _compressibility_via_density.diff(_molar_density).subs(
    _molar_density, 0)

assert expr_equals(_second_virial_coefficient, law.rhs)


@validate_input(
    bonding_forces_parameter_=attractive_forces_parameter,
    molecules_volume_parameter_=excluded_volume_parameter,
    temperature_=temperature,
)
@validate_output(second_virial_coefficient)
def calculate_second_virial_coefficient(
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
    temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        attractive_forces_parameter: bonding_forces_parameter_,
        excluded_volume_parameter: molecules_volume_parameter_,
        temperature: temperature_,
    })
    return Quantity(result)
