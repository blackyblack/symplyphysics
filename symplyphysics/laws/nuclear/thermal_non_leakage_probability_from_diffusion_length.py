from sympy import (Eq, solve, S)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    dimensionless,
    convert_to,
    validate_input,
)
from symplyphysics.core.symbols.probability import Probability

# Description
## Ptnl (thermal non-leakage factor) is the ratio of the number of thermal neutrons that do not leak from the
## reactor core during the neutron diffusion process to the number of neutrons that reach thermal energies.

## Law: Ptnl ≈ 1 / (1 + Lth^2 * Bg^2)
## Where:
## Lth - diffusion length of thermal neutrons.
##   See [diffusion area](./diffusion_area_from_diffusion_coefficient.py) implementation.
## Bg^2 - geometric buckling.
##   See [geometric buckling](./buckling/geometric_buckling_from_neutron_flux.py)
## Ptnl - thermal non-leakage probability.

thermal_diffusion_area = Symbol("thermal_diffusion_area", units.area)
geometric_buckling = Symbol("geometric_buckling", 1 / units.area)
thermal_non_leakage_probability = Symbol("thermal_non_leakage_probability", dimensionless)

law = Eq(thermal_non_leakage_probability, 1 / (1 + thermal_diffusion_area * geometric_buckling))


def print_law() -> str:
    return print_expression(law)


@validate_input(thermal_diffusion_area_=thermal_diffusion_area,
    geometric_buckling_=geometric_buckling)
def calculate_probability(thermal_diffusion_area_: Quantity,
    geometric_buckling_: Quantity) -> Probability:
    result_probability_expr = solve(law, thermal_non_leakage_probability,
        dict=True)[0][thermal_non_leakage_probability]
    result_expr = result_probability_expr.subs({
        thermal_diffusion_area: thermal_diffusion_area_,
        geometric_buckling: geometric_buckling_
    })
    result_factor = Quantity(result_expr)
    return Probability(convert_to(result_factor, S.One).evalf())
