from sympy import (Eq, Expr, solve, symbols, S, simplify)
from sympy.vector import Laplacian
from symplyphysics import (
    SI,
    Function,
    units,
    expr_to_quantity,
    Quantity,
    Symbol,
    print_expression,
    Dimensionless,
    convert_to,
    validate_input,
)
from symplyphysics.core.symbols.quantities import collect_factor_and_dimension

# Description
## The diffusion equation, based on Fick's law, provides an analytical solution of spatial neutron flux
## distribution in the multiplying system.

## Law: -D * Δ(Ф(x)) + Σa * Ф(x) = (1 / k) * v * Σf * Ф(x)
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

# NOTE: Neutron flux can be understood as scalar field. It is a function of the point in 3D space that returns a
# scalar value.

diffusion_coefficient = Symbol("diffusion_coefficient", units.length)
macroscopic_absorption_cross_section = Symbol("macroscopic_absorption_cross_section",
    1 / units.length)
macroscopic_fission_cross_section = Symbol("macroscopic_fission_cross_section", 1 / units.length)
effective_multiplication_factor = Symbol("effective_multiplication_factor", Dimensionless)
neutrons_per_fission = Symbol("neutrons_per_fission", Dimensionless)

# Position is a free variable of a function - do not specify its dimension
flux_position = symbols("flux_position")
neutron_flux = Function("neutron_flux", 1 / units.length**2 / units.time)
neutron_flux_laplacian = Function("neutron_flux_laplacian", 1 / units.length**4 / units.time)

neutron_flux_laplacian_definition = Eq(neutron_flux_laplacian(flux_position),
    Laplacian(neutron_flux(flux_position)),
    evaluate=False)

law = Eq(
    -1 * diffusion_coefficient * neutron_flux_laplacian(flux_position) +
    macroscopic_absorption_cross_section * neutron_flux(flux_position),
    (1 / effective_multiplication_factor) * neutrons_per_fission *
    macroscopic_fission_cross_section * neutron_flux(flux_position))

# As Laplacian is a second derivative over space coordinates (x, y, z), resulting dimension should be
# original dimension / units.length**2
assert neutron_flux_laplacian.dimension == neutron_flux.dimension / units.length**2


def print() -> str:
    return print_expression(law)


# neutron_flux_function_ should be a function on CoordSys3D
def apply_neutron_flux_function(neutron_flux_function_: Expr) -> Expr:
    # Manually divide to unit_length to get Laplacian dimension. CoordSys3D coordinates are dimensionless, hence
    # Laplacian cannot properly calculate resulting dimension.
    unit_length = Quantity(1, dimension=units.length)
    neutron_flux_laplacian_eval = neutron_flux_laplacian_definition.rhs.subs(
        neutron_flux(flux_position), neutron_flux_function_).doit() / unit_length**2
    applied_law = law.subs(neutron_flux_laplacian(flux_position), neutron_flux_laplacian_eval)
    applied_law = applied_law.subs(neutron_flux(flux_position), neutron_flux_function_)
    return simplify(applied_law)


# neutron_flux_function_ should be a function on CoordSys3D
# neutron_flux_function_ geometry should be defined with Quantity, eg width.dimension == units.length
@validate_input(neutrons_per_fission_=neutrons_per_fission,
    macroscopic_fission_cross_section_=macroscopic_fission_cross_section,
    macroscopic_absorption_cross_section_=macroscopic_absorption_cross_section,
    diffusion_coefficient_=diffusion_coefficient)
def calculate_multiplication_factor(neutron_flux_function_: Expr, neutrons_per_fission_: float,
    macroscopic_fission_cross_section_: Quantity, macroscopic_absorption_cross_section_: Quantity,
    diffusion_coefficient_: Quantity) -> float:

    (_, neutron_flux_dimension) = collect_factor_and_dimension(neutron_flux_function_)
    assert SI.get_dimension_system().equivalent_dims(neutron_flux_dimension, neutron_flux.dimension)

    applied_law = apply_neutron_flux_function(neutron_flux_function_)
    result_expr = applied_law.subs({
        neutrons_per_fission: neutrons_per_fission_,
        macroscopic_fission_cross_section: macroscopic_fission_cross_section_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_,
        diffusion_coefficient: diffusion_coefficient_
    })
    result_factor_expr = solve(result_expr, effective_multiplication_factor,
        dict=True)[0][effective_multiplication_factor]
    result_factor = expr_to_quantity(result_factor_expr)
    return convert_to(result_factor, S.One).evalf()
