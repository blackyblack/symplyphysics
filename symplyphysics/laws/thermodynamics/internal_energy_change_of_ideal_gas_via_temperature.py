"""
Internal energy change of ideal gas via temperature
===================================================

The internal energy of an ideal gas depends solely on its temperature and the number of gas particles
and is independent of other thermodynamic quantities such as pressure or density.

**Notation:**

#. :math:`d` denotes an exact, path-independent differential.
"""

from sympy import Eq
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

internal_energy_change = Symbol("internal_energy_change", units.energy)
"""
Infinitesimal change in internal energy of the system.

Symbol:
    :code:`dU`
"""

isochoric_heat_capacity = Symbol("isochoric_heat_capacity", units.energy / units.temperature)
r"""
Heat capacity at constant volume.

Symbol:
    :code:`C_V`

Latex:
    :math:`C_V`
"""

temperature_change = clone_symbol(symbols.thermodynamics.temperature, "temperature_change")
"""
Infinitesimal change in temperature of the system.

Symbol:
    :code:`dT`
"""

law = Eq(internal_energy_change, isochoric_heat_capacity * temperature_change)
r"""
:code:`dU = C_V * dT`

Latex:
    .. math::
        dU = C_V dT
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
