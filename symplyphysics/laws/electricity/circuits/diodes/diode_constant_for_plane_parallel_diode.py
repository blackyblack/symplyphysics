"""
Diode constant for parallel-plane vacuum diode
==============================================

The current-voltage characteristic of a vacuum diode is described by the
:math:`3/2`-power law (see :ref:`Current from voltage and diode constant in vacuum
diode`). The diode constant in this law depends only on the relative position, shape and
size of the electrodes of the vacuum diode. The vacuum diode may have different
electrode geometries. Here we are talking about a vacuum diode with parallel-plane
electrodes.

**Notation:**

#. :quantity_notation:`vacuum_permittivity`.
#. :quantity_notation:`elementary_charge`.
#. :quantity_notation:`electron_rest_mass`.

**Links:**

#. `Wikipedia, derivable from first formula <https://en.wikipedia.org/wiki/Space_charge#In_vacuum_(Child's_law)>`__.
"""

from sympy import Eq, Rational, solve, sqrt
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.quantities import vacuum_permittivity, elementary_charge, electron_rest_mass

diode_constant = symbols.diode_constant
"""
:symbols:`diode_constant`.
"""

electrode_area = symbols.area
"""
:symbols:`area` of an electrode.
"""

electrode_distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between electrodes.
"""

law = Eq(
    diode_constant,
    Rational(4, 9) * vacuum_permittivity * sqrt(2 * elementary_charge / electron_rest_mass) *
    (electrode_area / electrode_distance**2))


@validate_input(electrode_area_=electrode_area,
    distance_between_electrodes_=electrode_distance)
@validate_output(diode_constant)
def calculate_diode_constant(electrode_area_: Quantity,
    distance_between_electrodes_: Quantity) -> Quantity:
    result_expr = solve(law, diode_constant, dict=True)[0][diode_constant]
    result_expr = result_expr.subs({
        electrode_area: electrode_area_,
        electrode_distance: distance_between_electrodes_,
    })
    return Quantity(result_expr)
