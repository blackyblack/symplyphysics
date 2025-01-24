"""
Quality factor via energy loss
==============================

:math:`Q` factor (quality factor) is a dimensionless parameter that describes how
underdamped an oscillator or resonator is: the larger the :math:`Q` factor is, the less
damped it is. There are two nearly equivalent definitions of it that become
approximately equivalent as :math:`Q` becomes larger, meaning that the resonator becomes
less damped.

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

resonant_angular_frequency = clone_as_symbol(symbols.angular_frequency, subscript="\\text{r}")
"""
Resonant :symbols:`angular_frequency` of the oscillator.
"""

energy_stored = clone_as_symbol(symbols.energy, subscript="\\text{stored}")
"""
:symbols:`energy` stored in the oscillator.
"""

power_loss = clone_as_symbol(symbols.power, subscript="\\text{loss}")
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
