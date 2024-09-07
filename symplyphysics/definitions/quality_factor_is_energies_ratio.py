"""
Quality factor is ratio of energies
===================================

*Quality factor* is a property of an oscillating system defined as the ratio between the amount of
energy stored in system and the power losses.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, SymbolNew, validate_input, validate_output, dimensionless,
    angle_type)

quality_factor = SymbolNew("Q", dimensionless)
"""
Quality factor of the oscillator.
"""

resonant_frequency = SymbolNew("w", angle_type / units.time, display_latex="\\omega")
"""
Resonant angular frequency of the oscillator.
"""

stored_energy = SymbolNew("E", units.energy)
"""
Energy stored in the oscillator.
"""

dissipated_power = SymbolNew("P", units.power)
"""
Power dissipated from the oscillater.
"""

definition = Eq(quality_factor, resonant_frequency * stored_energy / dissipated_power)
"""
:laws:symbol::

:laws:latex::
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
