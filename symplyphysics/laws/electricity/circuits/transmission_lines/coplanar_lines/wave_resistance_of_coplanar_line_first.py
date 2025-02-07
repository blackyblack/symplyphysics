"""
Wave impedance of coplanar line (first)
=======================================

The coplanar transmission line is a dielectric substrate on the surface of which 3
electrodes are located. The wave resistance of a transmission line is a value determined
by the ratio of the voltage of the incident wave to the current of this wave in the
transmission line. When a wave propagates along a coplanar line, part of the field goes
out, since the coplanar line does not have metal borders on all sides, unlike, for
example, rectangular waveguides.

**Notes:**

#. Imagine an environment in which the field will have the same magnitude as the field
   of a microstrip line. The permittivity of such a medium will be called the effective
   permittivity of the line.

**Conditions:**

See the first and third conditions of :ref:`Effective permittivity of coplanar
transmission line when distance is greater than thickness`.

..
    TODO: rename file to feature wave *impedance*
    TODO: find link
"""

from sympy import Eq, solve, sqrt, root,  pi, log, evaluate
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
Effective :symbols:`relative_permittivity` of the coplanar line.
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
