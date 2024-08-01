"""
Quality factor is ratio of energies
===================================

*Quality factor* is a property of an oscillating system defined as the ratio between the amount of
energy stored in system and the power losses.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, dimensionless, angle_type)

quality_factor = Symbol("quality_factor", dimensionless)
"""
Quality factor of the oscillator.

Symbol:
    :code:`Q`
"""

resonant_frequency = Symbol("resonant_frequency", angle_type / units.time)
r"""
Resonant angular frequency of the oscillator.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

stored_energy = Symbol("stored_energy", units.energy)
"""
Energy stored in the oscillator.

Symbol:
    :code:`E`
"""

dissipated_power = Symbol("dissipated_power", units.power)
"""
Power dissipated from the oscillater.

Symbol:
    :code:`P`
"""

definition = Eq(quality_factor, resonant_frequency * stored_energy / dissipated_power)
r"""
:code:`Q = w * E / P`

Latex:
    .. math::
        Q = \frac{\omega E}{P}
"""


@validate_input(frequency_=resonant_frequency, energy_=stored_energy, power_=dissipated_power)
@validate_output(quality_factor)
def calculate_quality_factor(frequency_: Quantity, energy_: Quantity, power_: Quantity) -> Quantity:
    result_factor_expr = solve(definition, quality_factor, dict=True)[0][quality_factor]
    result_expr = result_factor_expr.subs({
        resonant_frequency: frequency_,
        stored_energy: energy_,
        dissipated_power: power_
    })
    return Quantity(result_expr)
