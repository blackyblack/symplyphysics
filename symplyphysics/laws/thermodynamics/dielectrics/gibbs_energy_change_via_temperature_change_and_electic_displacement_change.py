"""
Gibbs energy change via temperature change and electric displacement change
===========================================================================

The infinitesimal change in Gibbs energy of a system with a dielectric medium
can be expressed using the change in temperature and the change in electric
field strength.

**Conditions:**

#. The dielectric is isotropic whether or not the electric field is present.
#. The medium is homogeneous.
#. The volume change of the medium is insignificant.

**Links:**

#. Formula 31.8 on p. 122 of "General Course of Physics" (Obschiy kurs fiziki), vol. 3 by Sivukhin D.V. (1979).
"""

from sympy import Eq
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    clone_as_symbol,
    validate_input,
    validate_output,
)

gibbs_energy_density_change = Symbol("dG", units.energy / units.volume)
"""
Infinitesimal change in :symbols:`gibbs_energy` of the system per unit :symbols:`volume`.
"""

entropy_density = Symbol("S", units.energy / units.temperature / units.volume)
"""
:symbols:`entropy` of the system per unit :symbols:`volume`.
"""

temperature_change = clone_as_symbol(
    symbols.temperature,
    display_symbol="dT",
    display_latex="dT",
)
"""
Infinitesimal change in :symbols:`temperature` of the system.
"""

electric_displacement = symbols.electric_displacement
"""
:symbols:`electric_displacement`.
"""

electric_field_change = clone_as_symbol(
    symbols.electric_field_strength,
    display_symbol="dE",
    display_latex="dE",
)
"""
Infinitesimal change in :symbols:`electric_field_strength`.
"""

law = Eq(
    gibbs_energy_density_change,
    -1 * entropy_density * temperature_change - electric_displacement * electric_field_change,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derivation
## Make use of the formula for Gibbs energy and that for free energy.


@validate_input(
    entropy_density_=entropy_density,
    temperature_change_=temperature_change,
    electric_displacement_=electric_displacement,
    electric_field_change_=electric_field_change,
)
@validate_output(gibbs_energy_density_change)
def calculate_gibbs_energy_density_change(
    entropy_density_: Quantity,
    temperature_change_: Quantity,
    electric_displacement_: Quantity,
    electric_field_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        entropy_density: entropy_density_,
        temperature_change: temperature_change_,
        electric_displacement: electric_displacement_,
        electric_field_change: electric_field_change_,
    })
    return Quantity(result)
