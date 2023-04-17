from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Migration area (M^2) is equal to one-sixth of the square of the average distance (in all dimensions) between
## the neutron’s birth point (as a fast neutron) and its absorption (as a thermal neutron).

## Law: M^2 = L^2 + τth
## Where:
## L^2 - diffusion area.
##   See [diffusion area](./diffusion_area_from_diffusion_coefficient.py) implementation.
## τth - neutron Fermi age.
##   The Fermi age is related to the distance traveled during moderation, just as the diffusion length is for
##   thermal neutrons. The Fermi age is the same quantity as the slowing-down length squared, Ls^2, but the
##   slowing-down length is the square root of the Fermi age, τth = Ls^2. 
## M^2 - migration area.

diffusion_area = Symbol("diffusion_area", units.length**2)
neutron_fermi_age = Symbol("neutron_fermi_age", units.length**2)
migration_area = Symbol("migration_area", units.length**2)

law = Eq(migration_area, diffusion_area + neutron_fermi_age)

def print(expr: Expr) -> str:
    symbols = [diffusion_area, neutron_fermi_age, migration_area]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(diffusion_area_=diffusion_area, neutron_fermi_age_=neutron_fermi_age)
@validate_output_symbol(migration_area)
def calculate_migration_area(diffusion_area_: Quantity, neutron_fermi_age_: Quantity) -> Quantity:
    result_area_expr = solve(law, migration_area, dict=True)[0][migration_area]
    result_expr = result_area_expr.subs({
        diffusion_area: diffusion_area_,
        neutron_fermi_age: neutron_fermi_age_})
    return expr_to_quantity(result_expr)
