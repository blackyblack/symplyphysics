"""
Gibbs energy via enthalpy
=========================

Gibbs energy is a thermodynamic potential that can be used to calculate the maximum amount of work
that may be performed by a thermodynamically closed system at constant temperature and pressure.
It provides a necessary condition for processes such as chemical reactions that may occur under these
conditions.

**Links:**

#. `Wikipedia, first formula <https://en.wikipedia.org/wiki/Gibbs_free_energy>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

gibbs_energy = symbols.gibbs_energy
"""
:symbols:`gibbs_energy` of the system.
"""

enthalpy = symbols.enthalpy
"""
:symbols:`enthalpy` of the system.
"""

entropy = symbols.entropy
"""
:symbols:`entropy` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

law = Eq(gibbs_energy, enthalpy - temperature * entropy)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(thermal_effect_=enthalpy, entropy_=entropy, temperature_=temperature)
@validate_output(gibbs_energy)
def calculate_isobaric_potential(thermal_effect_: Quantity, entropy_: Quantity,
    temperature_: Quantity) -> Quantity:
    result_expr = solve(law, gibbs_energy, dict=True)[0][gibbs_energy]
    result_expr = result_expr.subs({
        enthalpy: thermal_effect_,
        entropy: entropy_,
        temperature: temperature_
    })
    return Quantity(result_expr)
