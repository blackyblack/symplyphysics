"""
Geometric buckling from multiplication factors and diffusion area
=================================================================

Geometric buckling can be expressed using infinite and effective multiplication factors
and diffusion area of the neutrons.

**Links:**

#. `Wikipedia, first part of second equation <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#Material_Buckling>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_macroscopic_fission_cross_section_diffusion_coefficient as buckling_law
from symplyphysics.laws.nuclear import diffusion_area_from_diffusion_coefficient as diffusion_area_law
from symplyphysics.laws.nuclear import infinite_multiplication_factor_from_macroscopic_fission_cross_section as infinite_multiplication_factor_law

infinite_multiplication_factor = symbols.infinite_multiplication_factor
"""
:symbols:`infinite_multiplication_factor`.
"""

effective_multiplication_factor = symbols.effective_multiplication_factor
"""
:symbols:`effective_multiplication_factor`.
"""

diffusion_area = symbols.neutron_diffusion_area
"""
:symbols:`neutron_diffusion_area`.
"""

geometric_buckling = symbols.geometric_buckling
"""
:symbols:`geometric_buckling`.
"""

law = Eq(geometric_buckling,
    (infinite_multiplication_factor / effective_multiplication_factor - 1) / diffusion_area)
"""
:laws:symbol::

:laws:latex::
"""

## Derive the same law from the diffusion area law and another geometric buckling law

_buckling_eq1 = buckling_law.law.subs({
    buckling_law.geometric_buckling: geometric_buckling,
    buckling_law.effective_multiplication_factor: effective_multiplication_factor
})
_diffusion_area_eq2 = diffusion_area_law.law.subs({
    diffusion_area_law.diffusion_area:
        diffusion_area,
    diffusion_area_law.diffusion_coefficient:
    buckling_law.diffusion_coefficient,
    diffusion_area_law.macroscopic_absorption_cross_section:
    buckling_law.macroscopic_absorption_cross_section
})
_infinite_multiplication_factor_eq3 = infinite_multiplication_factor_law.law.subs({
    infinite_multiplication_factor_law.infinite_multiplication_factor:
        infinite_multiplication_factor,
    infinite_multiplication_factor_law.neutrons_per_fission:
    buckling_law.neutrons_per_fission,
    infinite_multiplication_factor_law.macroscopic_fission_cross_section:
    buckling_law.macroscopic_fission_cross_section,
    infinite_multiplication_factor_law.macroscopic_absorption_cross_section:
    buckling_law.macroscopic_absorption_cross_section
})

_derived_law = [_buckling_eq1, _diffusion_area_eq2, _infinite_multiplication_factor_eq3]

## Check the equivalence of 'law' and '_derived_law'
_derived_geometric_buckling_squared = solve(_derived_law, (geometric_buckling,
    buckling_law.diffusion_coefficient, buckling_law.macroscopic_fission_cross_section),
    dict=True)[0][geometric_buckling]
assert expr_equals(law.rhs, _derived_geometric_buckling_squared)


@validate_input(infinite_multiplication_factor_=infinite_multiplication_factor,
    effective_multiplication_factor_=effective_multiplication_factor,
    diffusion_area_=diffusion_area)
@validate_output(geometric_buckling)
def calculate_geometric_buckling_squared(infinite_multiplication_factor_: float,
    effective_multiplication_factor_: float, diffusion_area_: Quantity) -> Quantity:
    result_buckling_expr = solve(law, geometric_buckling,
        dict=True)[0][geometric_buckling]
    result_expr = result_buckling_expr.subs({
        infinite_multiplication_factor: infinite_multiplication_factor_,
        effective_multiplication_factor: effective_multiplication_factor_,
        diffusion_area: diffusion_area_
    })
    return Quantity(result_expr)
