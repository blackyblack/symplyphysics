"""
Helmholtz free energy via internal energy
=========================================

*Helmholtz free energy* is a thermodynamic potential that measures the useful work obtainable
from a closed thermodynamic system at a constant temperature (isothermal process).
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

helmholtz_free_energy = Symbol("helmholtz_free_energy", units.energy)
"""
Helmholtz free energy of the system.

Symbol:
    :code:`F`
"""

internal_energy = Symbol("internal_energy", units.energy)
"""
Internal energy of the system.

Symbol:
    :code:`U`
"""

temperature = symbols.temperature
"""
:attr:`~symplyphysics.symbols.temperature` of the system.
"""

entropy = Symbol("entropy", units.energy / units.temperature)
"""
Entropy of the system.

Symbol:
    :code:`S`
"""

law = Eq(helmholtz_free_energy, internal_energy - temperature * entropy)
r"""
:code:`F = U - T * S`

Latex:
    .. math::
        F = U - T S
"""


@validate_input(
    internal_energy_=internal_energy,
    temperature_=temperature,
    entropy_=entropy,
)
@validate_output(helmholtz_free_energy)
def calculate_helmholtz_free_energy(
    internal_energy_: Quantity,
    temperature_: Quantity,
    entropy_: Quantity,
) -> Quantity:
    # Note that internal energy and entropy are only known up to a constant

    result = law.rhs.subs({
        internal_energy: internal_energy_,
        temperature: temperature_,
        entropy: entropy_,
    })
    return Quantity(result)
