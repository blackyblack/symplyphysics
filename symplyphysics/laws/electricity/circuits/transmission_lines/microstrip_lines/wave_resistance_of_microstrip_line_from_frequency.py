"""
Surge impedance of microstrip line from frequency
=================================================

The frequency-dependent surge impedance of a microstrip line can be calculated from its
frequency-independent surge impedance and effective permittivity.

..
    TODO: rename file to feature *surge impedance*
    TODO: find link
"""

from sympy import Eq, solve, sqrt, evaluate
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

surge_impedance = symbols.surge_impedance
"""
:symbols:`surge_impedance` of the microstrip line when frequency dependence is taken
into account.
"""

independent_surge_impedance = clone_as_symbol(symbols.surge_impedance, display_symbol="Z_S0", display_latex="Z_{\\text{S}, 0}")
"""
:symbols:`surge_impedance` of the microstrip line when frequency dependence is omitted.
"""

effective_permittivity = clone_as_symbol(
    symbols.relative_permittivity,
    display_symbol="epsilon_eff",
    display_latex="\\varepsilon_\\text{eff}",
)
"""
:symbols:`relative_permittivity` of the microstrip line when frequency dependence is
taken into account. See :ref:`Effective permittivity of microstrip line`.
"""

independent_effective_permittivity = clone_as_symbol(
    symbols.relative_permittivity,
    display_symbol="epsilon_eff0",
    display_latex="\\varepsilon_{\\text{eff}, 0}",
)
"""
:symbols:`relative_permittivity` of the microstrip line when frequency dependence is
omitted. See :ref:`Effective permittivity of microstrip line`.
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _first_expression = sqrt(independent_effective_permittivity / effective_permittivity)
    _second_expression = (effective_permittivity - 1) / (independent_effective_permittivity - 1)

law = Eq(surge_impedance, independent_surge_impedance * _first_expression * _second_expression)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(wave_resistance_without_frequency_=independent_surge_impedance,
    effective_permittivity_=effective_permittivity,
    effective_permittivity_without_frequency_=independent_effective_permittivity)
@validate_output(surge_impedance)
def calculate_wave_resistance(wave_resistance_without_frequency_: Quantity,
    effective_permittivity_: float, effective_permittivity_without_frequency_: float) -> Quantity:
    result_expr = solve(law, surge_impedance, dict=True)[0][surge_impedance]
    result_expr = result_expr.subs({
        independent_surge_impedance: wave_resistance_without_frequency_,
        effective_permittivity: effective_permittivity_,
        independent_effective_permittivity: effective_permittivity_without_frequency_
    })
    return Quantity(result_expr)
