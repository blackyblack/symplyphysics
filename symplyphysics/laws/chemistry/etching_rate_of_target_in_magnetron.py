from sympy import Eq, solve
from sympy.physics.units import elementary_charge, avogadro_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
)

# Description
## The ions of the gas-discharge plasma in the magnetron fall on the target and knock the atoms out of it.
## The sputtering coefficient shows how many target atoms are knocked out of the target by a single ion.
## The etching rate is how many nanometers of the target substance are etched per unit of time. In other words,
## this is how much thinner a target becomes per unit of time.

## Law is: V =  J * M * Y / (q * po * Na), where
## V - target etching rate,
## J - density of the current of the ion flux incident on the target,
## M - molar mass of the target atom,
## Y - sputtering coefficient,
## q - elementary charge,
## po - target density,
## Na - avogadro constant.

etching_rate = Symbol("etching_rate", units.velocity)

ion_current_density = Symbol("ion_current_density", units.current / units.area)
molar_mass_of_target_atom = Symbol("molar_mass_of_target_atom",
    units.mass / units.amount_of_substance)
sputtering_coefficient = Symbol("sputtering_coefficient", dimensionless)
target_density = Symbol("target_density", units.mass / units.volume)

law = Eq(
    etching_rate, ion_current_density * molar_mass_of_target_atom * sputtering_coefficient /
    (elementary_charge * target_density * avogadro_constant))


@validate_input(ion_current_density_=ion_current_density,
    molar_mass_of_target_atom_=molar_mass_of_target_atom,
    sputtering_coefficient_=sputtering_coefficient,
    target_density_=target_density)
@validate_output(etching_rate)
def calculate_etching_rate(ion_current_density_: Quantity, molar_mass_of_target_atom_: Quantity,
    sputtering_coefficient_: float, target_density_: Quantity) -> Quantity:
    result_expr = solve(law, etching_rate, dict=True)[0][etching_rate]
    result_expr = result_expr.subs({
        ion_current_density: ion_current_density_,
        molar_mass_of_target_atom: molar_mass_of_target_atom_,
        sputtering_coefficient: sputtering_coefficient_,
        target_density: target_density_,
    })
    return Quantity(result_expr)
