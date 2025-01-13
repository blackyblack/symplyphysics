from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, dimensionless, validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_macroscopic_fission_cross_section_diffusion_coefficient as buckling_law
from symplyphysics.laws.nuclear import diffusion_area_from_diffusion_coefficient as diffusion_area_law
from symplyphysics.laws.nuclear import infinite_multiplication_factor_from_macroscopic_fission_cross_section as infinite_multiplication_factor_law

# Description
## The quantity Bg^2 is called the geometrical buckling of the reactor and depends only on the geometry.
## This geometric buckling law is derived from more generic law:
## See: [geometric buckling](./geometric_buckling_from_macroscopic_fission_cross_section_diffusion_coefficient.py) implementation.

## Law: Bg^2 = (k_infinite / k_effective - 1) / L^2
## Where:
## k_infinite - infinite multiplication factor.
##   See [infinite multiplication factor](./infinite_multiplication_factor.py)
## k_effective - effective multiplication factor.
##   See [effective multiplication factor](./effective_multiplication_factor.py)
## L^2 - diffusion area.
##   See [diffusion area](./diffusion_area_from_diffusion_coefficient.py) implementation.
## Bg^2 - geometric buckling.
##   See [geometric buckling](./geometric_buckling_from_neutron_flux.py) implementation.

# Links
## Wikipedia, first part of second equation <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#Material_Buckling>

infinite_multiplication_factor = Symbol("infinite_multiplication_factor", dimensionless)
effective_multiplication_factor = Symbol("effective_multiplication_factor", dimensionless)
diffusion_area = Symbol("diffusion_area", units.area)
geometric_buckling_squared = Symbol("geometric_buckling_squared", 1 / units.area)

law = Eq(geometric_buckling_squared,
    (infinite_multiplication_factor / effective_multiplication_factor - 1) / diffusion_area)

## Derive the same law from the diffusion area law and another geometric buckling law

buckling_eq1 = buckling_law.law.subs({
    buckling_law.geometric_buckling_squared: geometric_buckling_squared,
    buckling_law.effective_multiplication_factor: effective_multiplication_factor
})
diffusion_area_eq2 = diffusion_area_law.law.subs({
    diffusion_area_law.diffusion_area:
        diffusion_area,
    diffusion_area_law.diffusion_coefficient:
    buckling_law.diffusion_coefficient,
    diffusion_area_law.macroscopic_absorption_cross_section:
    buckling_law.macroscopic_absorption_cross_section
})
infinite_multiplication_factor_eq3 = infinite_multiplication_factor_law.law.subs({
    infinite_multiplication_factor_law.infinite_multiplication_factor:
        infinite_multiplication_factor,
    infinite_multiplication_factor_law.neutrons_per_fission:
    buckling_law.neutrons_per_fission,
    infinite_multiplication_factor_law.macroscopic_fission_cross_section:
    buckling_law.macroscopic_fission_cross_section,
    infinite_multiplication_factor_law.macroscopic_absorption_cross_section:
    buckling_law.macroscopic_absorption_cross_section
})

derived_law = [buckling_eq1, diffusion_area_eq2, infinite_multiplication_factor_eq3]

## Check the equivalence of 'law' and 'derived_law'
derived_geometric_buckling_squared = solve(derived_law, (geometric_buckling_squared,
    buckling_law.diffusion_coefficient, buckling_law.macroscopic_fission_cross_section),
    dict=True)[0][geometric_buckling_squared]
assert expr_equals(law.rhs, derived_geometric_buckling_squared)


@validate_input(infinite_multiplication_factor_=infinite_multiplication_factor,
    effective_multiplication_factor_=effective_multiplication_factor,
    diffusion_area_=diffusion_area)
@validate_output(geometric_buckling_squared)
def calculate_geometric_buckling_squared(infinite_multiplication_factor_: float,
    effective_multiplication_factor_: float, diffusion_area_: Quantity) -> Quantity:
    result_buckling_expr = solve(law, geometric_buckling_squared,
        dict=True)[0][geometric_buckling_squared]
    result_expr = result_buckling_expr.subs({
        infinite_multiplication_factor: infinite_multiplication_factor_,
        effective_multiplication_factor: effective_multiplication_factor_,
        diffusion_area: diffusion_area_
    })
    return Quantity(result_expr)
