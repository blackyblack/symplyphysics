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
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

entropy_change = clone_as_symbol(symbols.entropy, display_symbol="dS", display_latex="dS")
"""
Infinitesimal change in :symbols:`entropy` of system :math:`B`.
"""

heat_supplied_to_system = clone_as_symbol(symbols.heat, display_symbol="delta(Q)", display_latex="\\delta Q")
"""
Infinitesimal amount of :symbols:`heat` transferred to system :math:`B`.
"""

common_temperature = symbols.temperature
"""
Common :symbols:`temperature` of systems :math:`A` and :math:`B`.
"""

law = Eq(entropy_change, heat_supplied_to_system / common_temperature)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    infinitesimal_transfer_of_heat_=heat_supplied_to_system,
    common_temperature_=common_temperature,
)
@validate_output(entropy_change)
def calculate_infinitesimal_entropy_change(
    infinitesimal_transfer_of_heat_: Quantity,
    common_temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        heat_supplied_to_system: infinitesimal_transfer_of_heat_,
        common_temperature: common_temperature_,
    })
    return Quantity(result)
