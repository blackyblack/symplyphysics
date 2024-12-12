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

**Links:**

#. `Wikipedia, first formula <https://en.wikipedia.org/wiki/Electrical_resistance_and_conductance#Relation_to_resistivity_and_conductivity>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols)

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the conductor.
"""

resistivity = symbols.electrical_resistivity
"""
:symbols:`electrical_resistivity` of the material.
"""

length = symbols.length
"""
:symbols:`length` of the conductor.
"""

area = symbols.area
"""
Cross-sectional :symbols:`area` of the conductor.
"""

law = Eq(resistance, resistivity * length / area)
"""
:laws:symbol::

:laws:latex::
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
