r"""
Enthalpy is internal energy plus pressure energy
================================================

Enthalpy :math:`H` of a thermodynamic system is defined as the sum of its internal energy
:math:`U` and the product of its pressure :math:`p` and volume :math:`V`, which is sometimes
referred to as the pressure energy.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

enthalpy = symbols.enthalpy
"""
:symbols:`enthalpy` of the system.
"""

internal_energy = symbols.internal_energy
"""
:symbols:`internal_energy` of the system.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` inside the system.
"""

volume = symbols.volume
"""
:symbols:`volume` of the system.
"""

law = Eq(enthalpy, internal_energy + pressure * volume)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    internal_energy_=internal_energy,
    pressure_=pressure,
    volume_=volume,
)
@validate_output(enthalpy)
def calculate_enthalpy(
    internal_energy_: Quantity,
    pressure_: Quantity,
    volume_: Quantity,
) -> Quantity:
    # Note that technically the internal energy is only known up to a constant

    result = law.rhs.subs({
        internal_energy: internal_energy_,
        pressure: pressure_,
        volume: volume_,
    })
    return Quantity(result)
