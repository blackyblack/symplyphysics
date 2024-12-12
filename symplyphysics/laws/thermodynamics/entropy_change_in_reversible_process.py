r"""
Entropy change in reversible process
====================================

The second law of thermodynamics applied to a closed system and an idealized, reversible or quasistatic process,
states that in such a process of transfer of energy as heat to closed thermodynamic system :math:`B`, which
allows energy but not matter exchanges, from an auxiliary thermodynamic system :math:`A`, an infinitesimal
increment in the entropy of system :math:`B` is defined to result from an infinitesimal transfer of heat to
system :math:`B` divided by the common thermodynamic temperature of systems :math:`A` and :math:`B`.

**Notation:**

#. :math:`\delta` (:code:`delta`) denotes an inexact, path-dependent differential.
#. :math:`d` denotes an exact, path-independent differential.

**Notes:**

#. Also applicable to actually possible quasi-static irreversible processes without composition change

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Entropy#Reversible_process>`__.
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

entropy_change = Symbol("entropy_change", units.energy / units.temperature)
"""
Infinitesimal change in entropy of system :math:`B`.

Symbol:
    :code:`dS`
"""

heat = Symbol("heat", units.energy)
r"""
Infinitesimal amount of heat transferred to system :math:`B`.

Symbol:
    :code:`delta(Q)`

Latex:
    :math:`\delta Q`
"""

common_temperature = symbols.temperature
"""
Common :symbols:`temperature` of systems :math:`A` and :math:`B`.
"""

law = Eq(entropy_change, heat / common_temperature)
r"""
:code:`dS = delta(Q) / T`

Latex:
    .. math::
        dS = \frac{\delta Q}{T}
"""


@validate_input(
    infinitesimal_transfer_of_heat_=heat,
    common_temperature_=common_temperature,
)
@validate_output(entropy_change)
def calculate_infinitesimal_entropy_change(
    infinitesimal_transfer_of_heat_: Quantity,
    common_temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        heat: infinitesimal_transfer_of_heat_,
        common_temperature: common_temperature_,
    })
    return Quantity(result)
