"""
Quality factor via energy loss
==============================

The quality factor of the system can be defined using the energy stored in the system
and the power dissipated from the system. In electrical system, the loss of energy
occurs due to resistors; in mechanical systems, it occurs due to external forces.

**Notes:**

#. There is :ref:`another definition <Quality factor via bandwidth>` that is
   approximately equivalent to this one at high quality factor values.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Q_factor#Stored_energy_definition>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

quality_factor = symbols.quality_factor
"""
:symbols:`quality_factor` of the oscillator.
"""

resonant_angular_frequency = clone_as_symbol(symbols.angular_frequency, display_symbol="w_r", display_latex="\\omega_\\text{r}")
"""
Resonant :symbols:`angular_frequency` of the oscillator.
"""

energy_stored = clone_as_symbol(symbols.energy, display_symbol="E_stored", display_latex="E_\\text{stored}")
"""
:symbols:`energy` stored in the oscillator.
"""

power_loss = clone_as_symbol(symbols.power, display_symbol="P_loss", display_latex="P_\\text{loss}")
"""
:symbols:`power` lost by the oscillator, i.e. the rate of energy dissipating from it.
"""

law = Eq(quality_factor, resonant_angular_frequency * (energy_stored / power_loss))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    resonant_angular_frequency_=resonant_angular_frequency,
    energy_stored_=energy_stored,
    power_loss_=power_loss,
)
@validate_output(quality_factor)
def calculate_quality_factor(
    resonant_angular_frequency_: Quantity,
    energy_stored_: Quantity,
    power_loss_: Quantity,
) -> float:
    result = law.rhs.subs({
        resonant_angular_frequency: resonant_angular_frequency_,
        energy_stored: energy_stored_,
        power_loss: power_loss_,
    })
    return Quantity(result).scale_factor
