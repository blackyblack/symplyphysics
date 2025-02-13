"""
Attenuation coefficient in microstrip metal when thickness is less than width times :math:`2 \\pi`
==================================================================================================

Under the conditions described below, the attenuation coefficient of the microstrip
metal can be calculated from the surge impedance of the line, the surface resistance of
the metal, the effective permittivity of the substrate, and the physical dimensions of
the system.

**Conditions:**

#. :math:`h \\ge w_\\text{eff}`.
#. :math:`h < 2 \\pi w_\\text{eff}`.

Here, :math:`h` is substrate thickness, and :math:`w_\\text{eff}` is effective width of
the microstrip.

..
    TODO: find link
"""

from sympy import Eq, solve, pi, log, evaluate
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

attenuation_coefficient = symbols.attenuation_coefficient
"""
:symbols:`attenuation_coefficient` of the metal in the microstrip.
"""

surface_resistance = clone_as_symbol(symbols.electrical_resistance, display_symbol="R_s", display_latex="R_\\text{s}")
"""
:symbols:`electrical_resistance` of the surface of the microstrip metal.
"""

surge_impedance = symbols.surge_impedance
"""
:symbols:`surge_impedance` of the microstrip line.
"""

effective_width = clone_as_symbol(symbols.length, display_symbol="w_eff", display_latex="w_\\text{eff}")
"""
Effective width (see :symbols:`length`) of the microstrip line. See :ref:`Effective
width of microstrip line`.
"""

width = clone_as_symbol(symbols.length, display_symbol="w", display_latex="w")
"""
Width (see :symbols:`length`) of the microstrip line.
"""

thickness = clone_as_symbol(symbols.thickness, display_symbol="t", display_latex="t")
"""
:symbols:`thickness` of the microstrip line.
"""

substrate_thickness = symbols.thickness
"""
:symbols:`thickness` of the substrate.
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    expression_1 = (effective_width / substrate_thickness)**2
    expression_2 = 1.38 * surface_resistance / (substrate_thickness * surge_impedance)
    expression_3 = 1 + (substrate_thickness / effective_width) * (1 - (1.25 / pi) * (thickness / substrate_thickness) + (1.25 / pi) * log(2 * (substrate_thickness / thickness)))

law = Eq(attenuation_coefficient,
    expression_2 * ((32 - expression_1) / (32 + expression_1)) * expression_3)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(surface_resistance_=surface_resistance,
    wave_resistance_=surge_impedance,
    thickness_of_substrate_=substrate_thickness,
    effective_width_=effective_width,
    strip_thickness_=thickness,
    width_=width)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(surface_resistance_: Quantity, wave_resistance_: Quantity,
    thickness_of_substrate_: Quantity, effective_width_: Quantity, strip_thickness_: Quantity,
    width_: Quantity) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    if thickness_of_substrate_.scale_factor < effective_width_.scale_factor:
        raise ValueError(
            "The thickness of substrate must be greater than or equal to the effective width")
    if thickness_of_substrate_.scale_factor >= effective_width_.scale_factor * 2 * pi:
        raise ValueError(
            "The thickness of substrate must be less than the effective width * 2 * pi")

    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        surface_resistance: surface_resistance_,
        surge_impedance: wave_resistance_,
        substrate_thickness: thickness_of_substrate_,
        effective_width: effective_width_,
        thickness: strip_thickness_,
        width: width_,
    })
    return Quantity(result_expr)
