"""
Internal energy change of ideal gas via temperature
===================================================

The internal energy of an ideal gas depends solely on its temperature and the number of gas particles
and is independent of other thermodynamic quantities such as pressure or density.

**Notation:**

#. :math:`d` denotes an exact, path-independent differential.

**Links:**

#. `Wikipedia, equivalent formula <https://en.wikipedia.org/wiki/Internal_energy#Internal_energy_of_the_ideal_gas>`__.
"""

from sympy import Eq
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

internal_energy_change = clone_as_symbol(symbols.internal_energy, display_symbol="dU", display_latex="dU")
"""
Infinitesimal change in :symbols:`internal_energy` of the system.
"""

isochoric_heat_capacity = clone_as_symbol(symbols.heat_capacity, subscript="V")
"""
:symbols:`heat_capacity` at constant :symbols:`volume`.
"""

temperature_change = clone_as_symbol(symbols.temperature, display_symbol="dT", display_latex="dT")
"""
Infinitesimal change in :symbols:`temperature` of the system.
"""

law = Eq(internal_energy_change, isochoric_heat_capacity * temperature_change)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive law from Joule-Thompson effect


@validate_input(
    isochoric_heat_capacity_=isochoric_heat_capacity,
    temperature_change_=temperature_change,
)
@validate_output(internal_energy_change)
def calculate_internal_energy(
    isochoric_heat_capacity_: Quantity,
    temperature_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        isochoric_heat_capacity: isochoric_heat_capacity_,
        temperature_change: temperature_change_,
    })
    return Quantity(result)
