from sympy.core.singleton import S
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units, Probability,
    validate_input, expr_to_quantity, convert_to
)

# Description
## Ptnl (thermal non-leakage factor) is the ratio of the number of thermal neutrons that do not leak from the
## reactor core during the neutron diffusion process to the number of neutrons that reach thermal energies.

## Law: Ptnl ≈ 1 / (1 + Lth^2 * Bg^2)
## Where:
## Lth - diffusion length of thermal neutrons.
##   See [diffusion area](./diffusion_area_from_diffusion_coefficient.py) implementation.
## Bg^2 - geometric buckling.
## Ptnl - thermal non-leakage probability.

thermal_diffusion_area = symbols('thermal_diffusion_area')
geometric_buckling = symbols('geometric_buckling')
thermal_non_leakage_probability = symbols('thermal_non_leakage_probability')

law = Eq(thermal_non_leakage_probability, 1 / (1 + thermal_diffusion_area * geometric_buckling))

def print():
    return pretty(law, use_unicode=False)

@validate_input(thermal_diffusion_area_=units.length**2, geometric_buckling_=(1 / units.length**2))
def calculate_probability(thermal_diffusion_area_: Quantity, geometric_buckling_: Quantity) -> Probability:
    result_probability_expr = solve(law, thermal_non_leakage_probability, dict=True)[0][thermal_non_leakage_probability]
    result_expr = result_probability_expr.subs({
        thermal_diffusion_area: thermal_diffusion_area_,
        geometric_buckling: geometric_buckling_})
    result_factor = expr_to_quantity(result_expr, 'thermal_non_leakage_factor')
    return Probability(convert_to(result_factor, S.One).n())
