"""
Resonance escape probability from resonance absorption integral
===============================================================

**Conditions:**

#. The reactor is homogeneous.
#. There are weak fast absorptions.
#. The absorber is predominant.

**Links:**

#. `Wikipedia, article <https://en.wikipedia.org/wiki/Resonance_escape_probability#Effective_resonance_integral>`__.
#. `Wikipedia, third row in table <https://en.wikipedia.org/wiki/Six_factor_formula>`__.
"""

from sympy import Eq, solve, exp
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    dimensionless,
    convert_to_float,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.symbols.probability import Probability

absorber_number_density = symbols.number_density
"""
:symbols:`number_density` of atoms in the absorber.
"""

effective_resonance_integral = Symbol("J_eff", units.area, display_latex="J_\\text{eff}")
"""
Effective resonance integral characterizes the absorption of neutrons by a single
nucleus in the resonance region.
"""

lethargy_gain_per_scattering = Symbol("xi", dimensionless, display_latex="\\xi")
"""
Average lethargy gain per scattering event. Lethargy is defined as decrease in neutron
energy.
"""

moderator_macroscopic_scattering_cross_section = clone_as_symbol(
    symbols.macroscopic_cross_section,
    display_symbol="Sigma_s",
    display_latex="\\Sigma_\\text{s}",
)
"""
:symbols:`macroscopic_cross_section` of scattering in the moderator.
"""

resonance_escape_probability = symbols.resonance_escape_probability
"""
:symbols:`resonance_escape_probability`.
"""

law = Eq(
    resonance_escape_probability,
    exp(-1 * (absorber_number_density * effective_resonance_integral) /
    (lethargy_gain_per_scattering * moderator_macroscopic_scattering_cross_section)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(absorber_atomic_number_density_=absorber_number_density,
    effective_resonance_integral_=effective_resonance_integral,
    average_lethargy_change_=lethargy_gain_per_scattering,
    macroscopic_scattering_cross_section_moderator_=moderator_macroscopic_scattering_cross_section)
@validate_output(resonance_escape_probability)
def calculate_resonance_escape_probability(
        absorber_atomic_number_density_: Quantity, effective_resonance_integral_: Quantity,
        average_lethargy_change_: float,
        macroscopic_scattering_cross_section_moderator_: Quantity) -> Probability:

    result_factor_expr = solve(law, resonance_escape_probability,
        dict=True)[0][resonance_escape_probability]
    result_expr = result_factor_expr.subs({
        absorber_number_density:
            absorber_atomic_number_density_,
        effective_resonance_integral:
            effective_resonance_integral_,
        lethargy_gain_per_scattering:
            average_lethargy_change_,
        moderator_macroscopic_scattering_cross_section:
            macroscopic_scattering_cross_section_moderator_
    })
    return Probability(convert_to_float(result_expr))
