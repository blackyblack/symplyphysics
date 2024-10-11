"""
Quality factor is ratio of energies
===================================

**Quality factor** is a property of an oscillating system defined as the ratio between the amount of
energy stored in system and the power losses.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

quality_factor = symbols.quality_factor
"""
:symbols:`quality_factor` of the oscillator.
"""

resonant_angular_frequency = symbols.angular_frequency
"""
Resonant :symbols:`angular_frequency` of the oscillator.
"""

stored_energy = symbols.energy
"""
:symbols:`energy` stored in the oscillator.
"""

dissipated_power = symbols.power
"""
:symbols:`power` dissipated from the oscillator.
"""

definition = Eq(quality_factor, resonant_angular_frequency * stored_energy / dissipated_power)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(frequency_=resonant_angular_frequency,
    energy_=stored_energy,
    power_=dissipated_power)
@validate_output(quality_factor)
def calculate_quality_factor(frequency_: Quantity, energy_: Quantity, power_: Quantity) -> Quantity:
    result_factor_expr = solve(definition, quality_factor, dict=True)[0][quality_factor]
    result_expr = result_factor_expr.subs({
        resonant_angular_frequency: frequency_,
        stored_energy: energy_,
        dissipated_power: power_
    })
    return Quantity(result_expr)
