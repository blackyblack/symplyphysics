"""
Surge impedance of microstrip line when effective width is less than substrate thickness
========================================================================================

The microstrip line is a dielectric substrate on which a metal strip is applied. When a
wave propagates along a microstrip line, part of the field goes out, since the
microstrip line does not have metal borders on all sides, unlike, for example,
rectangular waveguides.

**Conditions:**

#. Effective width :math:`w_\\text{eff}` of the microstrip line should be less than or
   equal to thickness :math:`h` of the substrate.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt, log, evaluate
from symplyphysics import (
    Quantity,
    units,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

surge_impedance = symbols.surge_impedance
"""
:symbols:`surge_impedance` of the microstrip line.
"""

effective_permittivity = clone_as_symbol(
    symbols.relative_permittivity,
    display_symbol="epsilon_eff",
    display_latex="\\varepsilon_\\text{eff}",
)
"""
Effective :symbols:`relative_permittivity` of the microstrip line. See :ref:`Effective
permittivity of microstrip line`.
"""

substrate_thickness = symbols.thickness
"""
:symbols:`thickness` of the substrate.
"""

effective_width = clone_as_symbol(symbols.length, display_symbol="w_eff", display_latex="w_\\text{eff}")
"""
Effective width (see :symbols:`length`) of the microstrip line. See :ref:`Effective
width of microstrip line`.
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _expression = log(8 * substrate_thickness / effective_width + effective_width / (4 * substrate_thickness))
    _resistance_constant = Quantity(60 * units.ohm, display_symbol="Z_0")

law = Eq(surge_impedance, (_resistance_constant / sqrt(effective_permittivity)) * _expression)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(effective_permittivity_=effective_permittivity,
    thickness_of_substrate_=substrate_thickness,
    effective_width_=effective_width)
@validate_output(surge_impedance)
def calculate_resistance(effective_permittivity_: float, thickness_of_substrate_: Quantity,
    effective_width_: Quantity) -> Quantity:
    if thickness_of_substrate_.scale_factor < effective_width_.scale_factor:
        raise ValueError("The thickness of substrate must be not less than the effective width")
    result_expr = solve(law, surge_impedance, dict=True)[0][surge_impedance]
    result_expr = result_expr.subs({
        effective_permittivity: effective_permittivity_,
        substrate_thickness: thickness_of_substrate_,
        effective_width: effective_width_
    })
    return Quantity(result_expr)
