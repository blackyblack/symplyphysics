r"""
Boltzmann factor via state energy and temperature
=================================================

The *Boltzmann factor* is an exponential factor that appears in many formulas of statistical physics
and thermodynamics, e.g. the canonical partition function of a classical discrete system.

**Notation:**

#. :math:`\exp` is the exponential function.
"""

from sympy import Eq, exp
from symplyphysics import (
    quantities,
    units,
    dimensionless,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    convert_to_float,
    SymbolNew,
)

boltzmann_factor = SymbolNew("f", dimensionless)
"""
Boltzmann factor.
"""

energy_of_state = SymbolNew("E[i]", units.energy, display_latex="E_i")
"""
Energy of state :math:`i`.
"""

equilibrium_temperature = symbols.temperature
"""
Equilibrium :attr:`~symplyphysics.symbols.temperature` of the system.
"""

boltzmann_constant = quantities.boltzmann_constant
"""
:attr:`~symplyphysics.quantities.boltzmann_constant`
"""

definition = Eq(boltzmann_factor,
    exp(-1 * energy_of_state / (boltzmann_constant * equilibrium_temperature)))
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
