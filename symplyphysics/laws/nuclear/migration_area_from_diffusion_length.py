from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

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

diffusion_area = symbols('diffusion_area')
neutron_fermi_age = symbols('neutron_fermi_age')
migration_area = symbols('migration_area')

law = Eq(migration_area, diffusion_area + neutron_fermi_age)

def print():
    return pretty(law, use_unicode=False)

@validate_input(diffusion_area_=units.length**2, neutron_fermi_age_=units.length**2)
@validate_output(units.length**2)
def calculate_migration_area(diffusion_area_: Quantity, neutron_fermi_age_: Quantity) -> Quantity:
    result_area_expr = solve(law, migration_area, dict=True)[0][migration_area]
    result_expr = result_area_expr.subs({
        diffusion_area: diffusion_area_,
        neutron_fermi_age: neutron_fermi_age_})
    return expr_to_quantity(result_expr, 'neutron_migration_area')
