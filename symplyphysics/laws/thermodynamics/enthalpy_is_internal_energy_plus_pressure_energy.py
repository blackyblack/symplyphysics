r"""
Enthalpy is internal energy plus pressure energy
================================================

Enthalpy :math:`H` of a thermodynamic system is defined as the sum of its internal energy
:math:`U` and the product of its pressure :math:`p` and volume :math:`V`, which is sometimes
referred to as the pressure energy.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

enthalpy = Symbol("enthalpy", units.energy)
"""
Enthalpy of the system.

Symbol:
    :code:`H`
"""

internal_energy = Symbol("internal_energy", units.energy)
"""
Internal energy of the system.

Symbol:
    :code:`U`
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure inside the system.

Symbol:
    :code:`p`
"""

volume = Symbol("volume", units.volume)
"""
Volume of the system.

Symbol:
    :code:`V`
"""

law = Eq(enthalpy, internal_energy + pressure * volume)
r"""
:code:`H = U + p * V`

Latex
    .. math::
        H = U + p V
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
