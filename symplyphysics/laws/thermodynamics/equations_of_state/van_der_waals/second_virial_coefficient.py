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

# Description
## The second virial coefficient is a coefficient appearing in the virial equation of state of
## a gas which describes the deviation from ideal gas behaviour. This coefficient describes pair
## interaction between molecules of the substance and can be found if the pair interaction potential
## is known, but it can also be derived from the equation of state as a series expansion with respect
## to inverse molar volume or equivalently to molar density.

# Law: B = b - a / (R * T)
## B - second virial coefficient of the virial expansion
## a - parameter of the van der Waals equation of state representing the magnitude of intermolecular forces
## b - parameter of the van der Waals equation of state representing the effective molecular size
## R - molar gas constant
## T - absolute temperature

# Conditions
## - Gas density is small enough within the context of perturbation theory.

second_virial_coefficient = Symbol("second_virial_coefficient", units.volume / units.amount_of_substance)
bonding_forces_parameter = Symbol(
    "bonding_forces_parameter",
    units.pressure * (units.volume / units.amount_of_substance)**2,
)
molecules_volume_parameter = Symbol(
    "molecules_volume_parameter",
    units.volume / units.amount_of_substance,
)
temperature = symbols.thermodynamics.temperature

law = Eq(
    second_virial_coefficient,
    molecules_volume_parameter - bonding_forces_parameter / (units.molar_gas_constant * temperature)
)

# Derive from the van der Waals equation of state and the virial equation

_volume = compressibility_def.volume
_mole_count = compressibility_def.amount_of_substance
_molar_density = SymSymbol("molar_density")

_pressure_expr = solve(vdw_eqn.law, vdw_eqn.pressure)[0].subs({
    vdw_eqn.volume: _volume,
    vdw_eqn.temperature: temperature,
    vdw_eqn.amount_of_substance: _mole_count,
    vdw_eqn.bonding_forces_parameter: bonding_forces_parameter,
    vdw_eqn.molecules_volume_parameter: molecules_volume_parameter,
})

_compressibility_via_volume = compressibility_def.definition.rhs.subs({
    compressibility_def.pressure: _pressure_expr,
    compressibility_def.temperature: temperature,
})

_mole_count_expr = density_qty_law.law.rhs.subs({
    density_qty_law.volumetric_density: _molar_density,
    density_qty_law.volume: _volume,
})

_compressibility_via_density = _compressibility_via_volume.subs(_mole_count, _mole_count_expr).simplify()

# Using the formula for Taylor series coefficients
_second_virial_coefficient = _compressibility_via_density.diff(_molar_density, 1).subs(_molar_density, 0)

assert expr_equals(_second_virial_coefficient, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    bonding_forces_parameter_=bonding_forces_parameter,
    molecules_volume_parameter_=molecules_volume_parameter,
    temperature_=temperature,
)
@validate_output(second_virial_coefficient)
def calculate_second_virial_coefficient(
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
    temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        bonding_forces_parameter: bonding_forces_parameter_,
        molecules_volume_parameter: molecules_volume_parameter_,
        temperature: temperature_,
    })
    return Quantity(result)
