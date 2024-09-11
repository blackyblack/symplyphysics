"""
Gibbs energy via enthalpy
=========================

Gibbs energy is a thermodynamic potential that can be used to calculate the maximum amount of work
that may be performed by a thermodynamically closed system at constant temperature and pressure.
It provides a necessary condition for processes such as chemical reactions that may occur under these
conditions.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

gibbs_energy = Symbol("gibbs_energy", units.energy / units.amount_of_substance)
"""
Gibbs energy of the system.

Symbol:
    :code:`G`
"""

enthalpy = Symbol("enthalpy", units.energy / units.amount_of_substance)
"""
Enthalpy of the system.

Symbol:
    :code:`H`
"""

entropy = Symbol("entropy", units.energy / units.amount_of_substance / units.temperature)
"""
Entropy of the system.

Symbol:
    :code:`S`
"""

temperature = symbols.temperature
"""
:attr:`~symplyphysics.symbols.temperature` of the system.
"""

law = Eq(gibbs_energy, enthalpy - temperature * entropy)
r"""
:code:`G = H - T * S`

Latex:
    .. math::
        G = H - T S
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
