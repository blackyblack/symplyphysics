"""
Wave impedance of coplanar line when hyperbolic sine ratio squared is between :math:`0` and :math:`\\frac{1}{2}`
================================================================================================================

Under the conditions described below, the wave impedance of a coplanar line depends on
its effective permittivity and physical dimensions.

**Conditions:**

#. :math:`h < \\frac{d}{4}`
#. :math:`0 < \\left( \\frac{\\sinh{ \\left((\\pi l) / (4 h)\\right) }}{\\sinh{ \\left((\\pi d) / (4 h)\\right) }} \\right)^2 \\le \\frac{1}{2}`

See below for symbol descriptions.

..
    TODO: check if it is *wave impedance* or *surge (characteristic) impedance*
    TODO: rename file to feature wave *impedance*
    TODO: find link
"""

from sympy import Eq, solve, sqrt, root, log, evaluate
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

wave_impedance = symbols.wave_impedance
"""
:symbols:`wave_impedance` of the coplanar line.
"""

effective_permittivity = clone_as_symbol(symbols.relative_permittivity, display_symbol="epsilon_eff", display_latex="\\varepsilon_\\text{eff}")
"""
Effective :symbols:`relative_permittivity` of the coplanar line. See :ref:`Effective permittivity of coplanar line <effective_permittivity_coplanar_line_def>`.
"""

electrode_distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between the first and last electrodes.
"""

central_electrode_width = symbols.length
"""
Width (see :symbols:`length`) of the central electrode of the coplanar line.
"""

resistance_constant = Quantity(30 * units.ohm, display_symbol="R_0")
"""
Constant equal to :math:`30 \\, \\Omega` (:code:`30 Ohm`).
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _first_expression = resistance_constant / sqrt(effective_permittivity)
    _second_expression = root(1 - (central_electrode_width / electrode_distance)**2, 4)
    _third_expression = log(2 * (1 + _second_expression) / (1 - _second_expression))

law = Eq(wave_impedance, _first_expression * _third_expression)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(effective_permittivity_=effective_permittivity,
    distance_between_electrodes_=electrode_distance,
    central_electrode_width_=central_electrode_width)
@validate_output(wave_impedance)
def calculate_wave_resistance(effective_permittivity_: float,
    distance_between_electrodes_: Quantity, central_electrode_width_: Quantity) -> Quantity:
    if (central_electrode_width_.scale_factor / distance_between_electrodes_.scale_factor)**2 > 0.5:
        raise ValueError("k^2 must be less than or equal to the 0.5")

    result_expr = solve(law, wave_impedance, dict=True)[0][wave_impedance]
    result_expr = result_expr.subs({
        effective_permittivity: effective_permittivity_,
        electrode_distance: distance_between_electrodes_,
        central_electrode_width: central_electrode_width_
    })
    return Quantity(result_expr)
