"""
Electric flux through closed surface in linear material
=======================================================

The :ref:`Gauss's law <Electric flux through closed surface via total charge>` can be written for
linear materials, featuring free charges, i.e. such charges that are not bound to any nucleus and
can move freely, as opposed to bound charges that occur due to the material's polarization.

**Conditions:**

#. The material is linear, homogeneous, isotropic, and nondispersive.

**Links:**

#. `Wikipedia — Gauss's law <https://en.wikipedia.org/wiki/Gauss%27s_law#Equation_for_linear_materials>`__.
"""

from sympy import Eq
from symplyphysics import symbols, clone_as_symbol, Quantity, validate_input, validate_output

electric_flux = symbols.electric_flux
"""
:symbols:`electric_flux` through a closed surface :math:`S`.
"""

total_free_charge = clone_as_symbol(
    symbols.charge,
    display_symbol="q_free",
    display_latex="q_\\text{free}",
)
"""
Total free :symbols:`charge` contained in the volume bound by :math:`S`.
"""

absolute_permittivity = symbols.absolute_permittivity
"""
:symbols:`absolute_permittivity` of the medium in the volume bound by :math:`S`.
"""

law = Eq(electric_flux, total_free_charge / absolute_permittivity)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    total_free_charge_=total_free_charge,
    absolute_permittivity_=absolute_permittivity,
)
@validate_output(electric_flux)
def calculate_total_electric_flux(
    total_free_charge_: Quantity,
    absolute_permittivity_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        total_free_charge: total_free_charge_,
        absolute_permittivity: absolute_permittivity_,
    })

    return Quantity(result)


# UNIQUE_LAW_ID: 518
