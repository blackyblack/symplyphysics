r"""
Electromotive force induced in rotating rod
===========================================

Let a rod rotate in a uniform magnetic field. The plane of rotation is perpendicular
to the magnetic field lines. The axis of rotation passes through one of the ends of
the rod. A wire is connected at both ends of the rod so that it makes a contour.
Then the electromotive force induced at the ends of the rod depends on the magnitude
of the magnetic flux density, the rotation frequency and the length of the rod.

**Links:**

#. `Example 13.4.2 <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/13%3A_Electromagnetic_Induction/13.04%3A_Motional_Emf>`__.

**Conditions:**

#. The angular velocity of the rod is parallel to the magnetic field. This means that
   the rod is rotating in a plane perpendicular to the magnetic field.
#. The magnetic field is uniform.
#. The angular velocity of the rod is constant.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import (
    electromotive_force_induced_in_moving_contour as _emf_law,
    magnetic_flux_from_induction_and_area as _magnetic_flux_law,
)
from symplyphysics.definitions import (
    angular_speed_is_angular_distance_derivative as _angular_speed_law,
)

electromotive_force = Symbol("electromotive_force", units.voltage)
r"""
Electromotive force induced in the rod.

Symbol:
    :code:`E`

Latex:
    :math:`\mathcal{E}`
"""

magnetic_flux_density = Symbol("magnetic_flux_density", units.magnetic_density)
"""
Magnitude of magnetic flux density.

Symbol:
    :code:`B`
"""

angular_frequency = Symbol("angular_frequency", angle_type / units.time)
r"""
Angular frequency of rod's rotation.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

length = Symbol("length", units.length)
"""
Length of the rod.

Symbol:
    :code:`l`
"""

law = Eq(electromotive_force, (1 / 2) * magnetic_flux_density * angular_frequency * length**2)
r"""
:code:`E = 1/2 * B * w * l^2`

Latex:
    .. math::
        \mathcal{E} = \frac{1}{2} B \omega l^2
"""

# Derive law

_time = _emf_law.time

_angle = _angular_speed_law.angular_distance(_time)

# Area of the circular sector
_loop_area = length**2 * _angle / 2

_magnetic_flux = _magnetic_flux_law.law.rhs.subs({
    _magnetic_flux_law.induction: magnetic_flux_density,
    _magnetic_flux_law.area: _loop_area,
    _magnetic_flux_law.angle: 0,  # the normal to the contour is parallel to the magnetic flux density
})

_emf = _emf_law.law.rhs.subs({
    _emf_law.current_turn_count: 1,
    _emf_law.magnetic_flux(_time): _magnetic_flux,
}).doit().subs({
    _angle.diff(_time): angular_frequency,
})

# We're interested in the absolute values of the EMF
assert expr_equals(abs(_emf), abs(law.rhs))


@validate_input(magnetic_induction_=magnetic_flux_density,
    rotation_frequency_=angular_frequency,
    rod_length_=length)
@validate_output(electromotive_force)
def calculate_voltage(magnetic_induction_: Quantity, rotation_frequency_: Quantity,
    rod_length_: Quantity) -> Quantity:
    result_expr = solve(law, electromotive_force, dict=True)[0][electromotive_force]
    result_expr = result_expr.subs({
        magnetic_flux_density: magnetic_induction_,
        angular_frequency: rotation_frequency_,
        length: rod_length_
    })
    return Quantity(result_expr)
