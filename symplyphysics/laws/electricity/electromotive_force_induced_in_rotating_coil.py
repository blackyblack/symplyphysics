r"""
Electromotive force induced in rotating coil
============================================

Suppose a coil is being rotated around the axis that lies in the coil's cross section,
see `Figure`_ in a magnetic field under the conditions described below. Then an 
electromotive will be induced in the contour of the coil. Its amplitude
depends on the number of turns in the coil, the magnetic flux density, the
angular frequency of the coil's rotation and the area of the coil's contour.

.. _Figure: https://www.schoolphysics.co.uk/age16-19/Electricity%20and%20magnetism/Electromagnetic%20induction/text/Induced_emf_in_a_rotating_coil/index.html

**Notes:**

#. The angle :math:`\varphi` between the normal to the coil's contour and the magnetic flux
   density is :math:`\varphi \propto \cos(\omega t)`.

**Conditions:**

#. The magnetic field is uniform.
#. The angular velocity of the coil's rotation is orthogonal to the magnetic field.
#. The area of the coil's contour is constant.
#. The rotation of the coil is uniform.
"""

from sympy import (Eq, solve, sin)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, dimensionless, angle_type)

electromotive_force = Symbol("electromotive_force", units.voltage)
r"""
Electromotive force induced in the coil.

Symbol:
    :code:`E`

Latex:
    :math:`\mathcal{E}`
"""

coil_turn_count = Symbol("coil_turn_count", dimensionless)
"""
Number of turns in the coil.

Symbol:
    :code:`N`
"""

magnetic_flux_density = Symbol("magnetic_flux_density", units.magnetic_density)
"""
Magnetic flux density.

Symbol:
    :code:`B`
"""

contour_area = Symbol("contour_area", units.area)
"""
Cross-sectional area of the contour enclosed by the coil.

Symbol:
    :code:`A`
"""

angular_frequency = Symbol("angular_frequency", angle_type / units.time)
r"""
Angular frequency of the coil's rotation.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

time = Symbol("time", units.time)

law = Eq(
    electromotive_force,
    -1 * coil_turn_count * magnetic_flux_density * contour_area * angular_frequency * sin(angular_frequency * time))
r"""
:code:`E = -1 * N * B * A * w * sin(w * t)`

Latex:
    .. math::
        \mathcal{E} = -N B A \omega \sin(\omega t)
"""


@validate_input(number_turns_=coil_turn_count,
    induction_=magnetic_flux_density,
    contour_area_=contour_area,
    rotation_frequency_=angular_frequency,
    time_=time)
@validate_output(electromotive_force)
def calculate_voltage(number_turns_: int, induction_: Quantity, contour_area_: Quantity,
    rotation_frequency_: Quantity, time_: Quantity) -> Quantity:
    result_expr = solve(law, electromotive_force, dict=True)[0][electromotive_force]
    result_expr = result_expr.subs({
        coil_turn_count: number_turns_,
        magnetic_flux_density: induction_,
        contour_area: contour_area_,
        angular_frequency: rotation_frequency_,
        time: time_
    })
    return Quantity(result_expr)
