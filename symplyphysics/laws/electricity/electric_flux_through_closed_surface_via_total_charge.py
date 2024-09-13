"""
Electric flux through closed surface via total charge
=====================================================

Known as **Gauss's law**, it states that the electric flux through any closed surface :math:`S` is equal to
the total charge within volume :math:`V` enclosed by that surface, and in the SI units, divided by the
vacuum permittivity. The closed surface is also referred to as **Gaussian surface**.

**Conditions:**

#. :math:`S = \partial V`, i.e. surface :math:`S` encloses volume :math:`V`. In other words, the surface
   must be closed.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    quantities,
    symbols,
)

total_electric_flux = symbols.electric_flux
"""
Electric flux through surface :math:`S`.
"""

total_charge = symbols.charge
"""
Total charge inside volume :math:`V`.
"""

law = Eq(total_electric_flux, total_charge / quantities.vacuum_permittivity)
r"""
:laws:symbol::

:laws:latex::
"""

@validate_input(total_charge_=total_charge)
@validate_output(total_electric_flux)
def calculate_total_electric_flux(total_charge_: Quantity) -> Quantity:
    result = law.rhs.subs(total_charge, total_charge_)
    return Quantity(result)
