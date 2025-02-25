"""
Diffusion equation from neutron flux
====================================

The diffusion equation, based on Fick's law, provides an analytical solution of spatial
neutron flux distribution in the multiplying system.

**Links:**

#. `NuclearPower, possible similar formula <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/diffusion-equation/>`__.
"""

from sympy import Eq, Expr, solve, simplify
from sympy.vector import Laplacian
from symplyphysics import (
    SI,
    Function,
    units,
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_function,
    clone_as_symbol,
)
from symplyphysics.core.dimensions import collect_factor_and_dimension

diffusion_coefficient = symbols.neutron_diffusion_coefficient
"""
:symbols:`neutron_diffusion_coefficient`.
"""

macroscopic_fission_cross_section = clone_as_symbol(symbols.macroscopic_cross_section, display_symbol="Sigma_f", display_latex="\\Sigma_text{f}")
"""
:symbols:`macroscopic_cross_section` of fission.
"""

macroscopic_absorption_cross_section = clone_as_symbol(symbols.macroscopic_cross_section, display_symbol="Sigma_a", display_latex="\\Sigma_\\text{a}")
"""
:symbols:`macroscopic_cross_section` of absorption.
"""

effective_multiplication_factor = symbols.effective_multiplication_factor
"""
:symbols:`effective_multiplication_factor`.
"""

neutrons_per_fission = clone_as_symbol(symbols.particle_count, display_symbol="nu", display_latex="\\nu")

# Position is a free variable of a function - do not specify its dimension
position = symbols.position
"""
:symbols:`position`.
"""

neutron_flux = clone_as_function(symbols.neutron_flux, [position])
"""
:symbols:`neutron_flux` as a function of :attr:`~position`.
"""

neutron_flux_laplacian = Function("Laplace(Phi)", [position], 1 / (units.length**4 * units.time), display_latex="\\nabla^{2} \\Phi")
"""
Laplacian of the :attr:`~neutron_flux` as a function of :attr:`~position`.
"""

neutron_flux_laplacian_definition = Eq(neutron_flux_laplacian(position),
    Laplacian(neutron_flux(position)),
    evaluate=False)

law = Eq(
    -1 * diffusion_coefficient * neutron_flux_laplacian(position) +
    macroscopic_absorption_cross_section * neutron_flux(position),
    (1 / effective_multiplication_factor) * neutrons_per_fission *
    macroscopic_fission_cross_section * neutron_flux(position))
"""
:laws:symbol::

:laws:latex::
"""

# As Laplacian is a second derivative over space coordinates (x, y, z), resulting dimension should be
# original dimension / units.length**2
assert SI.get_dimension_system().equivalent_dims(neutron_flux_laplacian.dimension,
    neutron_flux.dimension / units.length**2)


# neutron_flux_function_ should be a function on CoordSys3D
def apply_neutron_flux_function(neutron_flux_function_: Expr) -> Expr:
    # Manually divide to unit_length to get Laplacian dimension. CoordSys3D coordinates are dimensionless, hence
    # Laplacian cannot properly calculate resulting dimension.
    unit_length = Quantity(1, dimension=units.length)
    neutron_flux_laplacian_eval = neutron_flux_laplacian_definition.rhs.subs(
        neutron_flux(position), neutron_flux_function_).doit() / unit_length**2
    applied_law = law.subs(neutron_flux_laplacian(position), neutron_flux_laplacian_eval)
    applied_law = applied_law.subs(neutron_flux(position), neutron_flux_function_)
    return simplify(applied_law)


# neutron_flux_function_ should be a function on CoordSys3D
# neutron_flux_function_ geometry should be defined with Quantity, eg width.dimension == units.length
@validate_input(neutrons_per_fission_=neutrons_per_fission,
    macroscopic_fission_cross_section_=macroscopic_fission_cross_section,
    macroscopic_absorption_cross_section_=macroscopic_absorption_cross_section,
    diffusion_coefficient_=diffusion_coefficient)
@validate_output(effective_multiplication_factor)
def calculate_multiplication_factor(neutron_flux_function_: Expr, neutrons_per_fission_: float,
    macroscopic_fission_cross_section_: Quantity, macroscopic_absorption_cross_section_: Quantity,
    diffusion_coefficient_: Quantity) -> float:

    # this is like validate_input but does not require no free symbols
    (_, dimension) = collect_factor_and_dimension(neutron_flux_function_)
    assert SI.get_dimension_system().equivalent_dims(dimension, neutron_flux.dimension)

    applied_law = apply_neutron_flux_function(neutron_flux_function_)
    result_expr = applied_law.subs({
        neutrons_per_fission: neutrons_per_fission_,
        macroscopic_fission_cross_section: macroscopic_fission_cross_section_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_,
        diffusion_coefficient: diffusion_coefficient_
    })
    result_factor_expr = solve(result_expr, effective_multiplication_factor,
        dict=True)[0][effective_multiplication_factor]
    return convert_to_float(result_factor_expr)
