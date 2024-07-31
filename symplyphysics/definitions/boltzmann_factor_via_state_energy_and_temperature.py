r"""
Boltzmann factor via state energy and temperature
=================================================

The *Boltzmann factor* is an exponential factor that appears in many formulas of statistical physics
and thermodynamics, e.g. the canonical partition function of a classical discrete system.

**Notation:**

#. :math:`\exp` is the exponential function.
#. :math:`k_\text{B}` (:code:`k_B`) is the Boltzmann constant.
"""

from sympy import Eq, exp
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
    convert_to_float,
)

boltzmann_factor = Symbol("boltzmann_factor", dimensionless)
"""
Boltzmann factor.

Symbol:
    :code:`f`
"""

energy_of_state = Symbol("energy_of_state", units.energy)
"""
Energy of state :math:`i`.

Symbol:
    :code:`E_i`

Latex:
    :math:`E_i`
"""

equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature,
    "equilibrium_temperature")
"""
Equilibrium :attr:`~symplyphysics.symbols.thermodynamics.temperature` of the system.

Symbol:
    :code:`T`
"""

definition = Eq(boltzmann_factor,
    exp(-1 * energy_of_state / (units.boltzmann_constant * equilibrium_temperature)))
r"""
:code:`f = exp(-1 * E_i / (k_B * T))`

Latex:
    .. math::
        f = \exp{\left( - \frac{E_i}{k_\text{B} T} \right)}
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
