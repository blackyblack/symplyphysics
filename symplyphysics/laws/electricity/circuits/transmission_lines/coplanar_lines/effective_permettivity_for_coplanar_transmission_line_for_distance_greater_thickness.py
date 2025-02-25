"""
Effective permittivity of coplanar transmission line when distance is greater than thickness
============================================================================================

Under the conditions described below, the effective permittivity of a coplanar line can
be calculated from the relative permittivity of the substrate and the physical
dimensions of the system.

**Conditions:**

#. :math:`h < \\frac{d}{4}`
#. :math:`0 < \\left( \\frac{l}{d} \\right)^2 \\le \\frac{1}{2}`
#. :math:`0 < \\left( \\frac{\\sinh{ \\left((\\pi l) / (4 h)\\right) }}{\\sinh{ \\left((\\pi d) / (4 h)\\right) }} \\right)^2 \\le \\frac{1}{2}`

See below for symbol descriptions.

..
    TODO: fix file name
    TODO: add link
"""

from sympy import Eq, solve, root, pi, log, sinh, evaluate
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

effective_permittivity = clone_as_symbol(symbols.relative_permittivity, display_symbol="epsilon_eff", display_latex="\\varepsilon_\\text{eff}")
"""
Effective :symbols:`relative_permittivity` of the coplanar line. See :ref:`Effective permittivity of coplanar line <effective_permittivity_coplanar_line_def>`.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the dielectric substrate of the coplanar line.
"""

electrode_distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between the first and last electrodes.
"""

substrate_thickness = symbols.thickness
"""
:symbols:`thickness` of the substrate.
"""

central_electrode_width = symbols.length
"""
Width (see :symbols:`length`) of the central electrode of the coplanar line.
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _first_expression = (relative_permittivity - 1) / 2
    _second_expression = (
        sinh(pi * central_electrode_width / (4 * substrate_thickness))
        / sinh(pi * electrode_distance / (4 * substrate_thickness))
    )
    _third_expression = root(1 - _second_expression**2, 4)
    _fourth_expression = log(2 * (1 + _third_expression) / (1 - _third_expression))
    _fifth_expression = root(1 - (central_electrode_width / electrode_distance)**2, 4)
    _sixth_expression = log(2 * (1 + _fifth_expression) / (1 - _fifth_expression))
    _seventh_expression = _sixth_expression / _fourth_expression

law = Eq(effective_permittivity, 1 + _first_expression * _seventh_expression)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    distance_between_electrodes_=electrode_distance,
    thickness_of_substrate_=substrate_thickness,
    central_electrode_width_=central_electrode_width)
@validate_output(effective_permittivity)
def calculate_effective_permittivity(relative_permittivity_: float,
    distance_between_electrodes_: Quantity, thickness_of_substrate_: Quantity,
    central_electrode_width_: Quantity) -> float:
    if thickness_of_substrate_.scale_factor >= distance_between_electrodes_.scale_factor / 4:
        raise ValueError(
            "The thickness of substrate must be less than the distance between electrodes divided in 4"
        )
    if (central_electrode_width_.scale_factor / distance_between_electrodes_.scale_factor)**2 > 0.5:
        raise ValueError("k^2 must be less than or equal to the 0.5")
    if (sinh(pi * central_electrode_width_.scale_factor / 4 / thickness_of_substrate_.scale_factor)
            / sinh(pi * distance_between_electrodes_.scale_factor / 4 /
        thickness_of_substrate_.scale_factor))**2 > 0.5:
        raise ValueError("k1^2 must be less than or equal to the 0.5")
    result_expr = solve(law, effective_permittivity, dict=True)[0][effective_permittivity]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        electrode_distance: distance_between_electrodes_,
        substrate_thickness: thickness_of_substrate_,
        central_electrode_width: central_electrode_width_
    })
    return convert_to_float(result_expr)
