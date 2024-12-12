"""
Concentration of intrinsic charge carriers
==========================================

In the absence of external influences (lighting, electric field, etc.), there is a
nonzero concentration of free charge carriers in the semiconductor.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. There are no external influences, such as lighting, electric field, etc.

**Links:**

#. `Intrinsic carrier concentration <https://www.universitywafer.com/intrinsic-carrier-concentration.html>`_.
"""

from sympy import Eq, solve, exp, sqrt
from symplyphysics import (
    symbols,
    quantities,
    Quantity,
    validate_input,
    validate_output,
    clone_as_symbol,
)

charge_carriers_concentration = symbols.number_density
"""
:symbols:`number_density` of intrinsic charge carriers.
"""

density_of_states_in_conduction_band = clone_as_symbol(
    symbols.density_of_states,
    display_symbol="N_c",
    display_latex="N_\\text{c}",
)
"""
Effective :symbols:`density_of_states` in the conduction band.
"""

density_of_states_in_valence_band = clone_as_symbol(
    symbols.density_of_states,
    display_symbol="N_v",
    display_latex="N_\\text{v}",
)
"""
Effective :symbols:`density_of_states` in the valence band.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the semiconductor.
"""

band_gap = symbols.band_gap
"""
:symbols:`band_gap` of the semiconductor.
"""

law = Eq(
    charge_carriers_concentration,
    sqrt(density_of_states_in_conduction_band * density_of_states_in_valence_band) * exp(-band_gap /
    (2 * quantities.boltzmann_constant * temperature)))
"""
:laws:symbol::

:laws:latex::
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
