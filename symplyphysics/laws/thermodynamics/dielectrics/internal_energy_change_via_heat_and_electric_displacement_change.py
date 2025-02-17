"""
Internal energy change via heat and electric displacement change
================================================================

Internal energy change of the system with a dielectric medium can be
expressed using the infinitesimal heat flowing in or out of the system
and the change in electric displacement.

**Conditions:**

#. The dielectric is isotropic whether or not the electric field is present.
#. The medium is homogeneous.
#. The volume change of the medium is insignificant.

**Links:**

#. Formula 31.2 on p. 122 of "General Course of Physics" (Obschiy kurs fiziki), vol. 3 by Sivukhin D.V. (1979).
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

internal_energy_density_change = Symbol("dU", units.energy / units.volume)
"""
Infinitesimal change in :symbols:`internal_energy` per unit :symbols:`volume`.
"""

heat_density = Symbol("delta(Q)", units.energy / units.volume, display_latex="\\delta Q")
"""
Small amount of :symbols:`heat` added to or taken from the system per unit :symbols:`volume`.
"""

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength` in the medium.
"""

electric_displacement_change = clone_as_symbol(
    symbols.electric_displacement,
    display_symbol="dD",
    display_latex="dD",
)
"""
Infinitesimal change in :symbols:`electric_displacement`.
"""

law = Eq(internal_energy_density_change,
    heat_density + electric_field_strength * electric_displacement_change)
"""
:laws:symbol::

:laws:latex::
"""

# Derivation
# 1. First law of thermodynamics: `delta(Q) = dU + delta(A)`
# 2. Work done by the dielectric is composed of two parts: `delta(A) = p*dV - V*dot(E, dD)`
#    We omit the former term and focus on the latter one (TODO add this law).
# 3. Divide both parts of the (1) by `V`.


@validate_input(
    heat_density_=heat_density,
    electric_field_strength_=electric_field_strength,
    electric_displacement_change_=electric_displacement_change,
)
@validate_output(internal_energy_density_change)
def calculate_internal_energy_density_change(
    heat_density_: Quantity,
    electric_field_strength_: Quantity,
    electric_displacement_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        heat_density: heat_density_,
        electric_field_strength: electric_field_strength_,
        electric_displacement_change: electric_displacement_change_,
    })
    return Quantity(result)
