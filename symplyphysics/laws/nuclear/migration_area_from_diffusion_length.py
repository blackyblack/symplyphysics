"""
Migration area from diffusion length and Fermi age
==================================================

Migration area of a neutron can be found as the sum of the diffusion area and the Fermi
age.

**Links:**

#. `NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/migration-length-migration-area/>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

diffusion_area = symbols.neutron_diffusion_area
"""
:symbols:`neutron_diffusion_area`.
"""

fermi_age = symbols.neutron_fermi_age
"""
:symbols:`neutron_fermi_age`.
"""

migration_area = symbols.migration_area
"""
:symbols:`migration_area`.
"""

law = Eq(migration_area, diffusion_area + fermi_age)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(diffusion_area_=diffusion_area, neutron_fermi_age_=fermi_age)
@validate_output(migration_area)
def calculate_migration_area(diffusion_area_: Quantity, neutron_fermi_age_: Quantity) -> Quantity:
    result_area_expr = solve(law, migration_area, dict=True)[0][migration_area]
    result_expr = result_area_expr.subs({
        diffusion_area: diffusion_area_,
        fermi_age: neutron_fermi_age_
    })
    return Quantity(result_expr)
