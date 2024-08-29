"""
Resistance via resistivity and dimensions
=========================================

In the ideal conditions described below, the resistivity of a conductor is proportional
to its length and the inverse of its cross-sectional area. The constant of proportionality
is called resistivity of the material. Unlike resistance, resistivity is an intrinsic
property of the material and does not depend on its geometry.

**Conditions:**

#. The cross section is uniform throughout the conductor.
#. The current flows uniformly.
#. The conductor is made of a single material.
#. The electric field and current density are parallel and constant everywhere.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

resistance = Symbol("resistance", units.impedance)
"""
Resistance of the conductor.

Symbol:
    :code:`R`
"""

resistivity = Symbol("resistivity", units.impedance * units.length)
r"""
Resistivity of the material.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

length = Symbol("length", units.length)
"""
Length of the conductor.

Symbol:
    :code:`l`
"""

area = Symbol("area", units.area)
"""
Cross-sectional area of the conductor.

Symbol:
    :code:`A`
"""

law = Eq(resistance, resistivity * length / area)
r"""
:code:`R = rho * l / A`

Latex:
    .. math::
        R = \rho \frac{l}{A}
"""


@validate_input(resistivity_=resistivity, wire_length_=length, cross_section_=area)
@validate_output(resistance)
def calculate_resistance(resistivity_: Quantity, wire_length_: Quantity,
    cross_section_: Quantity) -> Quantity:
    result_resistance_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_resistance_expr.subs({
        resistivity: resistivity_,
        length: wire_length_,
        area: cross_section_
    })
    return Quantity(result_expr)
