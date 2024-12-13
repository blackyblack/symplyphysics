from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
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

# Links:
## NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/migration-length-migration-area/>

diffusion_area = Symbol("diffusion_area", units.area)
neutron_fermi_age = Symbol("neutron_fermi_age", units.length**2)
migration_area = Symbol("migration_area", units.area)

law = Eq(migration_area, diffusion_area + neutron_fermi_age)


@validate_input(diffusion_area_=diffusion_area, neutron_fermi_age_=neutron_fermi_age)
@validate_output(migration_area)
def calculate_migration_area(diffusion_area_: Quantity, neutron_fermi_age_: Quantity) -> Quantity:
    result_area_expr = solve(law, migration_area, dict=True)[0][migration_area]
    result_expr = result_area_expr.subs({
        diffusion_area: diffusion_area_,
        neutron_fermi_age: neutron_fermi_age_
    })
    return Quantity(result_expr)
