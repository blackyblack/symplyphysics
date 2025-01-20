"""
Free energy change via temperature change and electric displacement change
==========================================================================

The infinitesimal change in (Helmholtz) free energy of a system with a dielectric
medium can be expressed using the change in temperature and the change in electric
displacement.

**Conditions:**

#. The dielectric is isotropic whether or not the electric field is present.
#. The medium is homogeneous.
#. The volume change of the medium is insignificant.

**Links:**

#. Formula 31.7 on p. 122 of "General Course of Physics" (Obschiy kurs fiziki), vol. 3 by Sivukhin D.V. (1979).
"""

from sympy import Eq
from symplyphysics import (
    symbols,
    units,
    Quantity,
    SymbolNew,
    clone_as_symbol,
    validate_input,
    validate_output,
)

free_energy_density_change = SymbolNew("dH", units.energy / units.volume)
"""
Infinitesimal change in :symbols:`helmholtz_free_energy` of the system
per unit :symbols:`volume`.
"""

entropy_density = SymbolNew("S", units.energy / units.temperature / units.volume)
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

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`.
"""

electric_displacement_change = clone_as_symbol(
    symbols.electric_displacement,
    display_symbol="dD",
    display_latex="dD",
)
"""
Infinitesimal change in :symbols:`electric_displacement` of the system.
"""

law = Eq(
    free_energy_density_change,
    -1 * entropy_density * temperature_change +
    electric_field_strength * electric_displacement_change,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derivation
## Make use of the formula for free energy and the first law of thermodynamics applied to dielectric media.


@validate_input(
    entropy_density_=entropy_density,
    temperature_change_=temperature_change,
    electric_field_strength_=electric_field_strength,
    electric_displacement_change_=electric_displacement_change,
)
@validate_output(free_energy_density_change)
def calculate_free_energy_density_change(
    entropy_density_: Quantity,
    temperature_change_: Quantity,
    electric_field_strength_: Quantity,
    electric_displacement_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        entropy_density: entropy_density_,
        temperature_change: temperature_change_,
        electric_field_strength: electric_field_strength_,
        electric_displacement_change: electric_displacement_change_,
    })
    return Quantity(result)
