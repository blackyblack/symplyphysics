from sympy import Expr
from sympy.functions import exp
from symplyphysics import (
    Eq, pretty, solve, units, S, Probability, expr_to_quantity, convert_to
)
from symplyphysics.core.quantity_decorator import validate_input_symbols
from symplyphysics.core.symbols.quantities import Dimensionless, Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

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

geometric_buckling = Symbol("geometric_buckling", 1 / units.length**2)
neutron_fermi_age = Symbol("neutron_fermi_age", units.length**2)
fast_non_leakage_probability = Symbol("fast_non_leakage_probability", Dimensionless)

law = Eq(fast_non_leakage_probability, exp(-1 * geometric_buckling * neutron_fermi_age))

def print(expr: Expr) -> str:
    symbols = [geometric_buckling, neutron_fermi_age, fast_non_leakage_probability]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(
        geometric_buckling_=geometric_buckling,
        neutron_fermi_age_=neutron_fermi_age)
def calculate_probability(geometric_buckling_: Quantity, neutron_fermi_age_: Quantity) -> Probability:
    result_probability_expr = solve(law, fast_non_leakage_probability, dict=True)[0][fast_non_leakage_probability]
    result_expr = result_probability_expr.subs({
        geometric_buckling: geometric_buckling_,
        neutron_fermi_age: neutron_fermi_age_})
    result_factor = expr_to_quantity(result_expr)
    return Probability(convert_to(result_factor, S.One).evalf())
