"""
Geometric buckling from macroscopic cross sections and diffusion coefficient
============================================================================

Geometrical buckling is a quantity describing the reactor which depends only on its
geometry. It can also be calculated from the macroscopic cross sections and the
diffusion coefficient as well as the effective multiplicative factor and the neutron
production rate.

**Links:**

#. `Wikipedia, second part of second equation <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#Material_Buckling>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.nuclear import diffusion_equation_from_neutron_flux as diffusion_equation_law
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux as buckling_law

neutrons_per_fission = clone_as_symbol(symbols.particle_count,
    display_symbol="nu",
    display_latex="\\nu")
"""
The average number of neutrons produced per fission. See :symbols:`particle_count`.
"""

effective_multiplication_factor = symbols.effective_multiplication_factor
"""
:symbols:`effective_multiplication_factor`.
"""

macroscopic_fission_cross_section = clone_as_symbol(symbols.macroscopic_cross_section,
    display_symbol="Sigma_f",
    display_latex="\\Sigma_\\text{f}")
"""
:symbols:`macroscopic_cross_section` of fission.
"""

macroscopic_absorption_cross_section = clone_as_symbol(symbols.macroscopic_cross_section,
    display_symbol="Sigma_a",
    display_latex="\\Sigma_\\text{a}")
"""
:symbols:`macroscopic_cross_section` of absorption.
"""

diffusion_coefficient = symbols.neutron_diffusion_coefficient
"""
:symbols:`neutron_diffusion_coefficient`.
"""

geometric_buckling = symbols.geometric_buckling
"""
:symbols:`geometric_buckling`.
"""

law = Eq(geometric_buckling,
    ((neutrons_per_fission / effective_multiplication_factor) * macroscopic_fission_cross_section -
    macroscopic_absorption_cross_section) / diffusion_coefficient)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from the diffusion equation and geometric buckling from neutron flux law

_diffusion_eq1 = diffusion_equation_law.law.subs({
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
_buckling_eq2 = buckling_law.law.subs({
    buckling_law.geometric_buckling: geometric_buckling,
    buckling_law.neutron_flux: diffusion_equation_law.neutron_flux,
    buckling_law.position: diffusion_equation_law.position,
    buckling_law.neutron_flux_laplacian: diffusion_equation_law.neutron_flux_laplacian
})

derived_law = [
    _diffusion_eq1,
    _buckling_eq2,
    diffusion_equation_law._neutron_flux_laplacian_definition  # pylint: disable=protected-access
]

## Check the equivalence of 'law' and 'derived_law'
_derived_geometric_buckling_squared = solve(derived_law,
    (geometric_buckling, diffusion_equation_law.neutron_flux(diffusion_equation_law.position),
    diffusion_equation_law.neutron_flux_laplacian(diffusion_equation_law.position)),
    dict=True)[0][geometric_buckling]
assert expr_equals(law.rhs, _derived_geometric_buckling_squared)


@validate_input(neutrons_per_fission_=neutrons_per_fission,
    effective_multiplication_factor_=effective_multiplication_factor,
    macroscopic_fission_cross_section_=macroscopic_fission_cross_section,
    macroscopic_absorption_cross_section_=macroscopic_absorption_cross_section,
    diffusion_coefficient_=diffusion_coefficient)
@validate_output(geometric_buckling)
def calculate_buckling(neutrons_per_fission_: float, effective_multiplication_factor_: float,
    macroscopic_fission_cross_section_: Quantity, macroscopic_absorption_cross_section_: Quantity,
    diffusion_coefficient_: Quantity) -> Quantity:
    result_buckling_expr = solve(law, geometric_buckling, dict=True)[0][geometric_buckling]
    result_expr = result_buckling_expr.subs({
        neutrons_per_fission: neutrons_per_fission_,
        effective_multiplication_factor: effective_multiplication_factor_,
        macroscopic_fission_cross_section: macroscopic_fission_cross_section_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_,
        diffusion_coefficient: diffusion_coefficient_
    })
    return Quantity(result_expr)
