r"""
Concentration of intrinsic charge carriers
==========================================

In the absence of external influences (lighting, electric field, etc.), there is a
nonzero concentration of free charge carriers in the semiconductor.

**Notation:**

#. :math:`k_\text{B}` (:code:`k_B`) is the Boltzmann constant.

**Links:**

#. `Intrinsic carrier concentration <https://www.universitywafer.com/intrinsic-carrier-concentration.html>`_.

**Conditions:**

#. There are no external influences, such as lighting, electric field, etc.
"""

from sympy import (Eq, solve, sqrt, exp)
from sympy.physics.units import boltzmann
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

charge_carriers_concentration = Symbol("charge_carriers_concentration", 1 / units.volume)
r"""
Concentration of intrinsic charge carriers.

Symbol:
    :code:`n`
"""

density_of_states_in_conduction_band = Symbol("density_of_states_in_conduction_band",
    1 / units.volume)
r"""
Effective density of states in the conduction band.

Symbol:
    :code:`N_c`

Latex:
    :math:`N_\text{c}`
"""

density_of_states_in_valence_band = Symbol("density_of_states_in_valence_band", 1 / units.volume)
r"""
Effective density of states in the valence band.

Symbol:
    :code:`N_v`

Latex:
    :math:`N_\text{v}`
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the semiconductor.
"""

band_gap = Symbol("band_gap", units.energy)
r"""
Band gap of the semiconductor.

Symbol:
    :code:`E_g`

Latex:
    :math:`E_\text{g}`
"""

law = Eq(
    charge_carriers_concentration,
    sqrt(density_of_states_in_conduction_band * density_of_states_in_valence_band) * exp(-band_gap /
    (2 * boltzmann * temperature)))
r"""
:code:`n = sqrt(N_c * N_v) * exp(-1 * E_g / (2 * k_B * T))`

Latex:
    .. math::
        n = \sqrt{N_\text{c} N_\text{v}} \exp \left( -\frac{E_\text{g}}{2 k_\text{B} T} \right)
"""


@validate_input(density_of_states_in_conduction_band_=density_of_states_in_conduction_band,
    density_of_states_in_valence_band_=density_of_states_in_valence_band,
    band_gap_=band_gap,
    temperature_=temperature)
@validate_output(charge_carriers_concentration)
def calculate_concentration(density_of_states_in_conduction_band_: Quantity,
    density_of_states_in_valence_band_: Quantity, band_gap_: Quantity,
    temperature_: Quantity) -> Quantity:
    result_expr = solve(law, charge_carriers_concentration,
        dict=True)[0][charge_carriers_concentration]
    result_expr = result_expr.subs({
        density_of_states_in_conduction_band: density_of_states_in_conduction_band_,
        density_of_states_in_valence_band: density_of_states_in_valence_band_,
        band_gap: band_gap_,
        temperature: temperature_,
    })
    return Quantity(result_expr)
