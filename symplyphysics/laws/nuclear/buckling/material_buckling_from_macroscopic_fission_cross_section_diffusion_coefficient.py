from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units, simplify,
    validate_input, validate_output, expr_to_quantity
)

from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_macroscopic_fission_cross_section_diffusion_coefficient as buckling_law

# Description
## Material buckling (Bm^2) describes the difference between neutron production and neutron absorption.
## Bm^2 equals to geometric buckling (Bg^2) for critical reactor (k_effective = 1).

## Law: Bm^2 = (v * Σf - Σa) / D
## Where:
## v - average number of neutrons produced per fission.
## Σf - macroscopic fission cross-section of the fuel.
##   See [macroscopic cross-section](./macroscopic_cross_section_from_free_mean_path.py) implementation.
## Σa - macroscopic absorption cross-section of the fuel.
## D - diffusion coefficient.
##   See [diffusion coefficient](./neutron_diffusion_coefficient_from_scattering_cross_section.py) implementation.
## Bm^2 - material buckling.

neutrons_per_fission = symbols('neutrons_per_fission')
macroscopic_fission_cross_section = symbols('macroscopic_fission_cross_section')
macroscopic_absorption_cross_section = symbols('macroscopic_absorption_cross_section')
diffusion_coefficient = symbols('diffusion_coefficient')
material_buckling_squared = symbols('material_buckling_squared')

law = Eq(material_buckling_squared,
    (neutrons_per_fission * macroscopic_fission_cross_section - macroscopic_absorption_cross_section) /
    diffusion_coefficient)

# Derive the same law from the geometric buckling and critical reactor condition

buckling_eq1 = buckling_law.law.subs({
    buckling_law.geometric_buckling_squared: material_buckling_squared,
    buckling_law.neutrons_per_fission: neutrons_per_fission,
    buckling_law.macroscopic_fission_cross_section: macroscopic_fission_cross_section,
    buckling_law.macroscopic_absorption_cross_section: macroscopic_absorption_cross_section,
    buckling_law.diffusion_coefficient: diffusion_coefficient
})
critical_condition_eq2 = Eq(buckling_law.effective_multiplication_factor, 1)

derived_law = [buckling_eq1, critical_condition_eq2]

## Check the equivalence of 'law' and 'derived_law'
derived_material_buckling_squared = solve(
    derived_law, (material_buckling_squared, buckling_law.effective_multiplication_factor),
    dict=True)[0][material_buckling_squared]
assert simplify(law.rhs) == simplify(derived_material_buckling_squared)

def print():
    return pretty(law, use_unicode=False)

@validate_input(
    macroscopic_fission_cross_section_=(1 / units.length),
    macroscopic_absorption_cross_section_=(1 / units.length),
    diffusion_coefficient_=units.length)
@validate_output(1 / units.length**2)
def calculate_buckling(
    neutrons_per_fission_: float,
    macroscopic_fission_cross_section_: Quantity,
    macroscopic_absorption_cross_section_: Quantity,
    diffusion_coefficient_: Quantity) -> Quantity:
    result_buckling_expr = solve(law, material_buckling_squared, dict=True)[0][material_buckling_squared]
    result_expr = result_buckling_expr.subs({
        neutrons_per_fission: neutrons_per_fission_,
        macroscopic_fission_cross_section: macroscopic_fission_cross_section_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_,
        diffusion_coefficient: diffusion_coefficient_})
    return expr_to_quantity(result_expr, 'neutron_material_buckling_squared')
