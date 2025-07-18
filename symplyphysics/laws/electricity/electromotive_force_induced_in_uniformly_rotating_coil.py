"""
Electromotive force induced in rotating coil
============================================

Suppose a coil is being rotated around the axis that lies in the coil's cross section
(see `Figure <https://www.schoolphysics.co.uk/age16-19/Electricity%20and%20magnetism/Electromagnetic%20induction/text/Induced_emf_in_a_rotating_coil/index.html>`__)
in a magnetic field under the conditions described below. Then an electromotive will
be induced in the contour of the coil. Its amplitude depends on the number of turns in
the coil, the magnetic flux density, the angular frequency of the coil's rotation and
the area of the coil's contour.

**Notes:**

#. The angle :math:`\\varphi` between the normal to the coil's contour and the magnetic flux
   density is :math:`\\varphi \\propto \\cos(\\omega t)`. See `Figure <https://thefactfactor.com/wp-content/uploads/2020/03/Self-Induction-17.png>`__. 

**Conditions:**

#. The magnetic field is uniform.
#. The angular velocity of the coil's rotation is orthogonal to the magnetic field.
#. The area of the coil's contour is constant.
#. The angular speed of the coil's rotation constant.
"""

from sympy import Eq, solve, sin
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import (
    electromotive_force_induced_in_moving_contour as _emf_law,
    magnetic_flux_from_induction_and_area as _magnetic_flux_law,
)
from symplyphysics.laws.kinematics import (
    angular_position_via_constant_angular_speed_and_time as _angle_law,)

electromotive_force = symbols.electromotive_force
"""
:symbols:`electromotive_force` induced in the coil.
"""

coil_turn_count = symbols.positive_number
"""
Number of turns in the coil. See :symbols:`positive_number`.
"""

magnetic_flux_density = symbols.magnetic_flux_density
"""
:symbols:`magnetic_flux_density`.
"""

contour_area = symbols.area
"""
Cross-sectional :symbols:`area` of the contour enclosed by the coil.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the coil's rotation.
"""

time = symbols.time
"""
:symbols:`time`.
"""

law = Eq(
    electromotive_force, -1 * coil_turn_count * magnetic_flux_density * contour_area *
    angular_frequency * sin(angular_frequency * time))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law

_time = _emf_law.time

_angle = _angle_law.law.rhs.subs({
    _angle_law.initial_angular_position: 0,
    _angle_law.angular_speed: angular_frequency,
    _angle_law.time: _time,
})

_magnetic_flux = _magnetic_flux_law.law.rhs.subs({
    _magnetic_flux_law.magnetic_flux_density: magnetic_flux_density,
    _magnetic_flux_law.area: contour_area,
    _magnetic_flux_law.angle: _angle,
})

_emf = _emf_law.law.rhs.subs({
    _emf_law.current_turn_count: coil_turn_count,
    _emf_law.magnetic_flux(_time): _magnetic_flux,
}).doit().subs(_emf_law.time, time)

# We're interested in the absolute values of the EMF
assert expr_equals(abs(_emf), abs(law.rhs))


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
