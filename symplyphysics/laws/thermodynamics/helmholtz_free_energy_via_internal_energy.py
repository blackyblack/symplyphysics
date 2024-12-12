"""
Helmholtz free energy via internal energy
=========================================

*Helmholtz free energy* is a thermodynamic potential that measures the useful work obtainable
from a closed thermodynamic system at a constant temperature (isothermal process).

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Helmholtz_free_energy#Definition>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

helmholtz_free_energy = symbols.helmholtz_free_energy
"""
:symbols:`helmholtz_free_energy` of the system.
"""

internal_energy = symbols.internal_energy
"""
:symbols:`internal_energy` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

entropy = symbols.entropy
"""
:symbols:`entropy` of the system.
"""

law = Eq(helmholtz_free_energy, internal_energy - temperature * entropy)
"""
:laws:symbol::

:laws:latex::
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
