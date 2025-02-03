"""
Material buckling from material cross sections and diffusion coefficient
========================================================================

Material buckling can be calculated from material fission and absorption cross sections
as well as diffusion coefficient and average number of neutrons produced per fission.

**Links:**

#. `Wikipedia, last formula in paragraph <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#Material_Buckling>`__.
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
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_macroscopic_fission_cross_section_diffusion_coefficient as buckling_law

neutrons_per_fission = clone_as_symbol(symbols.particle_count, display_symbol="nu", display_latex="\\nu")
"""
The average number of neutrons produced per fission. See :symbols:`particle_count`.
"""

macroscopic_fission_cross_section = clone_as_symbol(symbols.macroscopic_cross_section, display_symbol="Sigma_f", display_latex="\\Sigma_\\text{f}")
"""
:symbols:`macroscopic_cross_section` of fission.
"""

macroscopic_absorption_cross_section = clone_as_symbol(symbols.macroscopic_cross_section, display_symbol="Sigma_a", display_latex="\\Sigma_\\text{a}")
"""
:symbols:`macroscopic_cross_section` of absorption.
"""

diffusion_coefficient = symbols.neutron_diffusion_coefficient
"""
:symbols:`neutron_diffusion_coefficient`.
"""

material_buckling = symbols.material_buckling
"""
:symbols:`material_buckling`.
"""

law = Eq(material_buckling, (neutrons_per_fission * macroscopic_fission_cross_section -
    macroscopic_absorption_cross_section) / diffusion_coefficient)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from the geometric buckling and critical reactor condition

_buckling_eq1 = buckling_law.law.subs({
    buckling_law.geometric_buckling: material_buckling,
    buckling_law.neutrons_per_fission: neutrons_per_fission,
    buckling_law.macroscopic_fission_cross_section: macroscopic_fission_cross_section,
    buckling_law.macroscopic_absorption_cross_section: macroscopic_absorption_cross_section,
    buckling_law.diffusion_coefficient: diffusion_coefficient
})
_critical_condition_eq2 = Eq(buckling_law.effective_multiplication_factor, 1)

_derived_law = [_buckling_eq1, _critical_condition_eq2]

## Check the equivalence of 'law' and '_derived_law'
_derived_material_buckling_squared = solve(_derived_law,
    (material_buckling, buckling_law.effective_multiplication_factor),
    dict=True)[0][material_buckling]
assert expr_equals(law.rhs, _derived_material_buckling_squared)


@validate_input(neutrons_per_fission_=neutrons_per_fission,
    macroscopic_fission_cross_section_=macroscopic_fission_cross_section,
    macroscopic_absorption_cross_section_=macroscopic_absorption_cross_section,
    diffusion_coefficient_=diffusion_coefficient)
@validate_output(material_buckling)
def calculate_buckling(neutrons_per_fission_: float, macroscopic_fission_cross_section_: Quantity,
    macroscopic_absorption_cross_section_: Quantity, diffusion_coefficient_: Quantity) -> Quantity:
    result_buckling_expr = solve(law, material_buckling,
        dict=True)[0][material_buckling]
    result_expr = result_buckling_expr.subs({
        neutrons_per_fission: neutrons_per_fission_,
        macroscopic_fission_cross_section: macroscopic_fission_cross_section_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_,
        diffusion_coefficient: diffusion_coefficient_
    })
    return Quantity(result_expr)
