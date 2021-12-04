from sympy import simplify
from sympy.core.expr import Expr
from sympy.core.singleton import S
from sympy.vector import Laplacian
from symplyphysics import (
    symbols, Eq, solve, pretty, Function, Quantity, units,
    validate_input, expr_to_quantity, convert_to
)

# Description
## The diffusion equation, based on Fick's law, provides an analytical solution of spatial neutron flux
## distribution in the multiplying system.

## Law: -D * Δ(Ф(x)) + Σa * Ф(x) = (1 / k) * Σf * Ф(x)
## Where:
## D - diffusion coefficient.
##   See [diffusion coefficient](./neutron_diffusion_coefficient_from_scattering_cross_section.py) implementation.
## Δ - Laplacian operator.
## Ф(x) (neutron flux density) - number of neutrons crossing through some arbitrary cross-sectional unit area in all
##   directions per unit time.
## Σf - macroscopic fission cross-section of the fuel.
##   See [macroscopic cross-section](./macroscopic_cross_section_from_free_mean_path.py) implementation.
## Σa - macroscopic absorption cross-section of the fuel.
## k - effective multiplication factor.
##   See [effective multiplication factor](./effective_multiplication_factor.py)
## v - average number of neutrons produced per fission.
## x - coordinates for the neutron flux density.

diffusion_coefficient = symbols('diffusion_coefficient')
macroscopic_absorption_cross_section = symbols('macroscopic_absorption_cross_section')
macroscopic_fission_cross_section = symbols('macroscopic_fission_cross_section')
neutron_flux_function = symbols('neutron_flux_function', cls = Function)
effective_multiplication_factor = symbols('effective_multiplication_factor')
neutrons_per_fission = symbols('neutrons_per_fission')
flux_position = symbols('flux_position')

law = Eq(
    -1 * diffusion_coefficient * Laplacian(neutron_flux_function(flux_position)) +
    macroscopic_absorption_cross_section * neutron_flux_function(flux_position),
    (1 / effective_multiplication_factor) * neutrons_per_fission * macroscopic_fission_cross_section *
    neutron_flux_function(flux_position))

# rearrange the law and check it is the same as original
rearranged_law = Eq(
    -1 * Laplacian(neutron_flux_function(flux_position)) / neutron_flux_function(flux_position),
    ((1 / effective_multiplication_factor) * neutrons_per_fission * macroscopic_fission_cross_section - macroscopic_absorption_cross_section) /
    diffusion_coefficient)

assert solve(law, neutron_flux_function(flux_position))[0] == solve(rearranged_law, neutron_flux_function(flux_position))[0]

def print():
    return pretty(law, use_unicode=False)

# neutron_flux_function_ should be a function on CoordSys3D
def apply_neutron_flux_function(neutron_flux_function_: Function) -> Expr:
    applied_law = law.subs(neutron_flux_function(flux_position), neutron_flux_function_)
    return simplify(applied_law.doit())

# neutron_flux_function_ should be a function on CoordSys3D
# neutron_flux_function_ geometry should be defined with Quantity, eg width.dimension == units.length
@validate_input(
    macroscopic_fission_cross_section_=(1 / units.length),
    macroscopic_absorption_cross_section_=(1 / units.length),
    diffusion_coefficient_=units.length)
def calculate_multiplication_factor(
    neutron_flux_function_: Function,
    neutrons_per_fission_: float,
    macroscopic_fission_cross_section_: Quantity,
    macroscopic_absorption_cross_section_: Quantity,
    diffusion_coefficient_: Quantity) -> float:

    result_expr = law.subs({
        neutron_flux_function(flux_position): neutron_flux_function_,
        neutrons_per_fission: neutrons_per_fission_,
        macroscopic_fission_cross_section: macroscopic_fission_cross_section_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_,
        diffusion_coefficient: diffusion_coefficient_})

    result_factor_expr = solve(result_expr, effective_multiplication_factor, dict=True)[0][effective_multiplication_factor]
    result_factor = expr_to_quantity(result_factor_expr, 'eff_multiplication_factor')
    return convert_to(result_factor, S.One).n()
