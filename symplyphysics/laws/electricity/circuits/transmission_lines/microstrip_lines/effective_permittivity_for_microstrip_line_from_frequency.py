"""
Effective permittivity of microstrip line from frequency
========================================================

The frequency-dependent effective permittivity of the microstrip line can be calculated
from its frequency-indendent effective permittivity and physical dimensions.

**Notation:**

#. :quantity_notation:`speed_of_light`.

..
    TODO: find link
"""

from sympy import Eq, Rational, solve, sqrt, log, evaluate
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
    quantities,
)

effective_permittivity = clone_as_symbol(
    symbols.relative_permittivity,
    display_symbol="epsilon_eff",
    display_latex="\\varepsilon_\\text{eff}",
)
"""
Effective :symbols:`relative_permittivity` of the microstrip line when frequency
dependence is taken into account. See :ref:`Effective permittivity of microstrip line`.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the dielectric substrate of the microstrip line.
"""

frequency = symbols.temporal_frequency
"""
:symbols:`temporal_frequency` of the signal.
"""

substrate_thickness = symbols.thickness
"""
:symbols:`thickness` of the substrate.
"""

width = clone_as_symbol(symbols.length, display_symbol="w", display_latex="w")
"""
Width (see :symbols:`length`) of the microstrip line.
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
    _first_expression = (4 / (quantities.speed_of_light * 2)) * substrate_thickness * sqrt(relative_permittivity - 1)
    _second_expression = (1 + 2 * log(1 + width / substrate_thickness))**2
    _third_expression = _first_expression * _second_expression
    _fourth_expression = sqrt(relative_permittivity) - sqrt(independent_effective_permittivity)
    _fifth_expression = _fourth_expression / (1 + 4 * (_third_expression * frequency)**Rational(-3, 2))

law = Eq(
    effective_permittivity,
    (_fifth_expression + sqrt(independent_effective_permittivity))**2,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    frequency_=frequency,
    thickness_of_substrate_=substrate_thickness,
    width_=width,
    effective_permittivity_without_frequency_=independent_effective_permittivity)
@validate_output(effective_permittivity)
def calculate_effective_permittivity(relative_permittivity_: float, frequency_: Quantity,
    thickness_of_substrate_: Quantity, width_: Quantity,
    effective_permittivity_without_frequency_: float) -> float:
    result_expr = solve(law, effective_permittivity, dict=True)[0][effective_permittivity]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        frequency: frequency_,
        substrate_thickness: thickness_of_substrate_,
        width: width_,
        independent_effective_permittivity: effective_permittivity_without_frequency_
    })
    return convert_to_float(result_expr)
