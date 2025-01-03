"""
Enthalpy change via entropy change and electric field change
============================================================

The infinitesimal change in entropy of the system with a dielectric medium can
be expressed using the change in its entropy and the change in the electric
field strength.

**Conditions:**

#. The dielectric is isotropic whether or not the electric field is present.
#. The medium is homogeneous.
#. The volume change of the medium is insignificant.

**Links:**

#. Formula 31.9 on p. 122 of "General Course of Physics" (Obschiy kurs fiziki), vol. 3 by Sivukhin D.V. (1979).
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

enthalpy_density_change = SymbolNew("dH", units.energy / units.volume)
"""
Infinitesimal change in :symbols:`enthalpy` of the system per unit :symbols:`volume`.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

entropy_density_change = SymbolNew("dS", units.energy / units.temperature / units.volume)
"""
Infinitesimal change in :symbols:`entropy` of the system per unit :symbols:`volume`.
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
    enthalpy_density_change,
    temperature * entropy_density_change - electric_displacement * electric_field_change,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derivation
## Make use of the formula for enthalpy and that of internal energy.


@validate_input(
    temperature_=temperature,
    entropy_density_change_=entropy_density_change,
    electric_displacement_=electric_displacement,
    electric_field_change_=electric_field_change,
)
@validate_output(enthalpy_density_change)
def calculate_enthalpy_density_change(
    temperature_: Quantity,
    entropy_density_change_: Quantity,
    electric_displacement_: Quantity,
    electric_field_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        temperature: temperature_,
        entropy_density_change: entropy_density_change_,
        electric_displacement: electric_displacement_,
        electric_field_change: electric_field_change_,
    })
    return Quantity(result)
