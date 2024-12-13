from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, dimensionless, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.nuclear import diffusion_equation_from_neutron_flux as diffusion_equation_law
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux as buckling_law

# Description
## The quantity Bg^2 is called the geometrical buckling of the reactor and depends only on the geometry.
## See [geometric buckling](./geometric_buckling_from_neutron_flux.py) implementation.

## Law: Bg^2 = ((v / k_effective) * Σf - Σa) / D
## Where:
## v - average number of neutrons produced per fission.
## k_effective - effective multiplication factor.
##   See [effective multiplication factor](./effective_multiplication_factor.py)
## Σf - macroscopic fission cross-section of the fuel.
##   See [macroscopic cross-section](./macroscopic_cross_section_from_free_mean_path.py) implementation.
## Σa - macroscopic absorption cross-section of the fuel.
## D - diffusion coefficient.
##   See [diffusion coefficient](./neutron_diffusion_coefficient_from_scattering_cross_section.py) implementation.
## Bg^2 - geometric buckling.

# Links
## Wikipedia, second part of second equation <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#Material_Buckling>

neutrons_per_fission = Symbol("neutrons_per_fission", dimensionless)
effective_multiplication_factor = Symbol("effective_multiplication_factor", dimensionless)
macroscopic_fission_cross_section = Symbol("macroscopic_fission_cross_section", 1 / units.length)
macroscopic_absorption_cross_section = Symbol("macroscopic_absorption_cross_section",
    1 / units.length)
diffusion_coefficient = Symbol("diffusion_coefficient", units.length)
geometric_buckling_squared = Symbol("geometric_buckling_squared", 1 / units.area)

law = Eq(geometric_buckling_squared,
    ((neutrons_per_fission / effective_multiplication_factor) * macroscopic_fission_cross_section -
    macroscopic_absorption_cross_section) / diffusion_coefficient)

# Derive the same law from the diffusion equation and geometric buckling from neutron flux law

diffusion_eq1 = diffusion_equation_law.law.subs({
    diffusion_equation_law.effective_multiplication_factor:
        effective_multiplication_factor,
    diffusion_equation_law.diffusion_coefficient:
        diffusion_coefficient,
    diffusion_equation_law.macroscopic_absorption_cross_section:
        macroscopic_absorption_cross_section,
    diffusion_equation_law.macroscopic_fission_cross_section:
        macroscopic_fission_cross_section,
    diffusion_equation_law.neutrons_per_fission:
        neutrons_per_fission
})
buckling_eq2 = buckling_law.law.subs({
    buckling_law.geometric_buckling_squared: geometric_buckling_squared,
    buckling_law.neutron_flux: diffusion_equation_law.neutron_flux,
    buckling_law.flux_position: diffusion_equation_law.flux_position,
    buckling_law.neutron_flux_laplacian: diffusion_equation_law.neutron_flux_laplacian
})

derived_law = [
    diffusion_eq1, buckling_eq2, diffusion_equation_law.neutron_flux_laplacian_definition
]

## Check the equivalence of 'law' and 'derived_law'
derived_geometric_buckling_squared = solve(derived_law, (geometric_buckling_squared,
    diffusion_equation_law.neutron_flux(diffusion_equation_law.flux_position),
    diffusion_equation_law.neutron_flux_laplacian(diffusion_equation_law.flux_position)),
    dict=True)[0][geometric_buckling_squared]
assert expr_equals(law.rhs, derived_geometric_buckling_squared)


def print_law() -> str:
    return print_expression(law)


@validate_input(neutrons_per_fission_=neutrons_per_fission,
    effective_multiplication_factor_=effective_multiplication_factor,
    macroscopic_fission_cross_section_=macroscopic_fission_cross_section,
    macroscopic_absorption_cross_section_=macroscopic_absorption_cross_section,
    diffusion_coefficient_=diffusion_coefficient)
@validate_output(geometric_buckling_squared)
def calculate_buckling(neutrons_per_fission_: float, effective_multiplication_factor_: float,
    macroscopic_fission_cross_section_: Quantity, macroscopic_absorption_cross_section_: Quantity,
    diffusion_coefficient_: Quantity) -> Quantity:
    result_buckling_expr = solve(law, geometric_buckling_squared,
        dict=True)[0][geometric_buckling_squared]
    result_expr = result_buckling_expr.subs({
        neutrons_per_fission: neutrons_per_fission_,
        effective_multiplication_factor: effective_multiplication_factor_,
        macroscopic_fission_cross_section: macroscopic_fission_cross_section_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_,
        diffusion_coefficient: diffusion_coefficient_
    })
    return Quantity(result_expr)
