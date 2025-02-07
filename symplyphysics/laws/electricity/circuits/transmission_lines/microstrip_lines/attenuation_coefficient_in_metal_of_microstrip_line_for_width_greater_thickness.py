"""
Attenuation coefficient in metal of microstrip line when width is greater than thickness
========================================================================================

The microstrip line is a dielectric substrate on which a metal strip is applied. When a
wave propagates along a microstrip line, part of the field goes out, since the
microstrip line does not have metal borders on all sides, unlike, for example,
rectangular waveguides. Then imagine an environment in which the field will have the
same magnitude as the field of a microstrip line. The permittivity of such a medium will
be called the effective permittivity of the line.

**Conditions:**

#. The thickness of the substrate of the microstrip line should be less than the
   effective width.

..
    TODO: find link
"""

from sympy import Eq, solve, pi, log, evaluate
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

attenuation_coefficient = symbols.attenuation_coefficient
"""
:symbols:`attenuation_coefficient` of the metal of the microstrip line.
"""

surface_resistance = clone_as_symbol(symbols.electrical_resistance, display_symbol="R_s", display_latex="R_\\text{s}")
"""
:symbols:`electrical_resistance` of the surface of the metal strip.
"""

surge_impedance = symbols.surge_impedance
"""
:symbols:`surge_impedance` of the microstrip line.
"""

substrate_thickness = symbols.thickness
"""
:symbols:`thickness` of the substrate.
"""

effective_width = clone_as_symbol(symbols.length, display_symbol="w_eff", display_latex="w_\\text{eff}")
"""
Effective width (see :symbols:`length`) of the microstrip line. It is the width of such
a flat capacitor, the electric field between the plates of which is equal to the
electric field in the dielectric of the substrate under the line strip.
"""

thickness = clone_as_symbol(symbols.thickness, display_symbol="t", display_latex="t")
"""
:symbols:`thickness` of the strip of the microstrip line.
"""

effective_permittivity = clone_as_symbol(symbols.relative_permittivity, display_symbol="epsilon_eff", display_latex="\\varepsilon_\\text{eff}")
"""
Effective :symbols:`relative_permittivity` of the microstrip line.
"""

constant = Quantity(6.1e-5 / units.ohm**2, display_symbol="a")
"""
Constant equal to :math:`6.1 \\cdot 10^{-5} \\, \\Omega^{-2}` (:code:`6.1e-5 Ohm^(-2)`).
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _relative_width = effective_width / substrate_thickness
    expression_2 = constant * surface_resistance * surge_impedance * effective_permittivity / substrate_thickness
    expression_3 = 1 + (1 / _relative_width) * (1 - (1.25 / pi) * (thickness / substrate_thickness) + (1.25 / pi) * log(2 * substrate_thickness / thickness))

law = Eq(
    attenuation_coefficient,
    expression_2 * (_relative_width + 0.667 * _relative_width / (_relative_width + 1.444)) * expression_3)
"""
:laws:symbol::

:laws:latex::

..
    NOTE: check if the numbers in the formula should actually be `2/3` and `13/9`
"""


@validate_input(surface_resistance_=surface_resistance,
    wave_resistance_=surge_impedance,
    thickness_of_substrate_=substrate_thickness,
    effective_width_=effective_width,
    strip_thickness_=thickness,
    effective_permittivity_=effective_permittivity)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(surface_resistance_: Quantity, wave_resistance_: Quantity,
    thickness_of_substrate_: Quantity, effective_width_: Quantity, strip_thickness_: Quantity,
    effective_permittivity_: float) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    if thickness_of_substrate_.scale_factor >= effective_width_.scale_factor:
        raise ValueError("The thickness of substrate must be less than the effective width")

    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        surface_resistance: surface_resistance_,
        surge_impedance: wave_resistance_,
        substrate_thickness: thickness_of_substrate_,
        effective_width: effective_width_,
        thickness: strip_thickness_,
        effective_permittivity: effective_permittivity_,
    })
    return Quantity(result_expr)
