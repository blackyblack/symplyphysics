from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Material buckling (Bm^2) describes the characteristics of the fuel material in an infinite medium.

## Law: Bm^2 = (v * Σf - Σa) / (k_effective * D)
## Where:
## v - average number of neutrons produced per fission.
## Σf - macroscopic fission cross-section of the fuel.
##   See [macroscopic transport cross-section](./macroscopic_transport_cross_section.py) implementation.
## Σa - macroscopic absorption cross-section of the fuel.
## k_effective - effective multiplication factor.
##   See [effective multiplication factor](./effective_multiplication_factor.py) implementation.
## D - diffusion coefficient.
##   See [diffusion coefficient](./neutron_diffusion_coefficient_from_scattering_cross_section.py) implementation.
## Bm^2 - material buckling.

neutrons_per_fission = symbols('neutrons_per_fission')
macroscopic_fission_cross_section = symbols('macroscopic_fission_cross_section')
macroscopic_absorption_cross_section = symbols('macroscopic_absorption_cross_section')
effective_multiplication_factor = symbols('effective_multiplication_factor')
diffusion_coefficient = symbols('diffusion_coefficient')
material_buckling = symbols('material_buckling')

law = Eq(material_buckling,
    (neutrons_per_fission * macroscopic_fission_cross_section - macroscopic_absorption_cross_section) /
    (effective_multiplication_factor * diffusion_coefficient))

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
    effective_multiplication_factor_: float,
    diffusion_coefficient_: Quantity) -> Quantity:
    result_buckling_expr = solve(law, material_buckling, dict=True)[0][material_buckling]
    result_expr = result_buckling_expr.subs({
        neutrons_per_fission: neutrons_per_fission_,
        macroscopic_fission_cross_section: macroscopic_fission_cross_section_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_,
        effective_multiplication_factor: effective_multiplication_factor_,
        diffusion_coefficient: diffusion_coefficient_})
    return expr_to_quantity(result_expr, 'neutron_material_buckling')
