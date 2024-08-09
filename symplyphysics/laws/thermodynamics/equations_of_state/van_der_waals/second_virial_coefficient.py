r"""
Second virial coefficient
=========================

The *second virial coefficient* is a coefficient appearing in the virial equation of state of
a gas which describes the deviation from ideal gas behaviour. This coefficient describes pair
interaction between molecules of the substance and can be found if the pair interaction potential
is known, but it can also be derived from the equation of state as a series expansion with respect
to inverse molar volume or equivalently, to molar density.

**Notation:**

#. :math:`R` is the molar gas constant.

**Conditions:**

#. The gas density is small enough within the context of perturbation theory.
"""

from sympy import Eq, solve, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import compressibility_factor_is_deviation_from_ideal_gas as compressibility_def
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import equation as vdw_eqn
from symplyphysics.laws.quantities import quantity_is_volumetric_density_times_volume as density_qty_law

second_virial_coefficient = Symbol("second_virial_coefficient",
    units.volume / units.amount_of_substance)
"""
Second virial coefficient of the virial expansion.

Symbol:
    :code:`B`
"""

attractive_forces_parameter = Symbol(
    "attractive_forces_parameter",
    units.pressure * (units.volume / units.amount_of_substance)**2,
)
"""
Parameter of the van der Waals equation denoting the magnitude of attractive
forces between gas molecules.

Symbol:
    :code:`a`
"""

excluded_volume_parameter = Symbol(
    "excluded_volume_parameter",
    units.volume / units.amount_of_substance,
)
"""
Parameter of the van der Waals equation denoting an excluded molar molar_volume
due to a finite size of molecules.

Symbol:
    :code:`b`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the van der Waals fluid.

Symbol:
    :code:`T`
"""

law = Eq(
    second_virial_coefficient, excluded_volume_parameter - attractive_forces_parameter /
    (units.molar_gas_constant * temperature))
r"""
:code:`B = b - a / (R * T)`

Latex:
    .. math::
        B = b - \frac{a}{R T}
"""

# Derive from the van der Waals equation of state and the virial equation

_volume = compressibility_def.volume
_mole_count = compressibility_def.amount_of_substance
_molar_density = SymSymbol("molar_density")

_pressure_expr = solve(vdw_eqn.law, vdw_eqn.pressure)[0].subs({
    vdw_eqn.molar_volume: _volume / _mole_count,
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


def print_law() -> str:
    return print_expression(law)


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
