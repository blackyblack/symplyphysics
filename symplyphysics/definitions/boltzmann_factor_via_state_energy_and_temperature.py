"""
Boltzmann factor via state energy and temperature
=================================================

The *Boltzmann factor* is an exponential factor that appears in many formulas of statistical physics
and thermodynamics, e.g. the canonical partition function of a classical discrete system.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Boltzmann_distribution#>`__.
"""

from sympy import Eq, exp
from symplyphysics import (
    quantities,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    convert_to_float,
    clone_as_symbol,
)

boltzmann_factor = symbols.boltzmann_factor
"""
:symbols:`boltzmann_factor`.
"""

energy_of_state = clone_as_symbol(symbols.energy, display_symbol="E[i]", display_latex="E_i")
"""
:symbols:`energy` of state :math:`i`.
"""

equilibrium_temperature = symbols.temperature
"""
Equilibrium :symbols:`temperature` of the system.
"""

definition = Eq(
    boltzmann_factor,
    exp(-1 * energy_of_state / (quantities.boltzmann_constant * equilibrium_temperature)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    energy_of_state_=energy_of_state,
    equilibrium_temperature_=equilibrium_temperature,
)
@validate_output(boltzmann_factor)
def calculate_boltzmann_factor(
    energy_of_state_: Quantity,
    equilibrium_temperature_: Quantity,
) -> float:
    result = definition.rhs.subs({
        energy_of_state: energy_of_state_,
        equilibrium_temperature: equilibrium_temperature_,
    })
    return convert_to_float(result)
