from sympy import (Eq, solve, exp, S)
from symplyphysics import (
    units,
    expr_to_quantity,
    Quantity,
    Symbol,
    print_expression,
    Dimensionless,
    convert_to,
    validate_input,
)
from symplyphysics.core.symbols.probability import Probability

# Description
## Ptnl (fast non-leakage factor) is the ratio of the number of fast neutrons that do not leak from the reactor
## core during the slowing down process to the number of fast neutrons produced by fissions at all energies.

## Law: Pfnl ≈ e^(-Bg^2 * τth)
## Where:
## e - exponent.
## Bg^2 - geometric buckling.
##   See [geometric buckling](./buckling/geometric_buckling_from_neutron_flux.py) implementation.
## τth - neutron Fermi age.
##   The Fermi age is related to the distance traveled during moderation, just as the diffusion length is for
##   thermal neutrons. The Fermi age is the same quantity as the slowing-down length squared (Ls^2).
## Pfnl - fast non-leakage probability.

geometric_buckling = Symbol("geometric_buckling", 1 / units.area)
neutron_fermi_age = Symbol("neutron_fermi_age", units.length**2)
fast_non_leakage_probability = Symbol("fast_non_leakage_probability", Dimensionless)

law = Eq(fast_non_leakage_probability, exp(-1 * geometric_buckling * neutron_fermi_age))


def print_law() -> str:
    return print_expression(law)


@validate_input(geometric_buckling_=geometric_buckling, neutron_fermi_age_=neutron_fermi_age)
def calculate_probability(geometric_buckling_: Quantity,
    neutron_fermi_age_: Quantity) -> Probability:
    result_probability_expr = solve(law, fast_non_leakage_probability,
        dict=True)[0][fast_non_leakage_probability]
    result_expr = result_probability_expr.subs({
        geometric_buckling: geometric_buckling_,
        neutron_fermi_age: neutron_fermi_age_
    })
    result_factor = expr_to_quantity(result_expr)
    return Probability(convert_to(result_factor, S.One).evalf())
