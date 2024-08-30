r"""
Electromotive force induced in rotating rod
===========================================

Let a rod rotate in a uniform magnetic field. The plane of rotation is perpendicular
to the magnetic field lines. The axis of rotation passes through one of the ends of
the rod. Then the electromotive force induced at the ends of the rod depends on the
magnitude of the magnetic induction, the rotation frequency and the length of the rod.

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

electromotive_force = Symbol("electromotive_force", units.voltage)
r"""
Electromotive force induced in the rod.

Symbol:
    :code:`E`

Latex:
    :math:`\mathcal{E}`
"""

magnetic_field = Symbol("magnetic_field", units.magnetic_flux_density)
"""
Value of the magnetic field.

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

law = Eq(electromotive_force, (1 / 2) * magnetic_field * angular_frequency * length**2)
r"""
:code:`E = 1/2 * B * w * l^2`

Latex:
    .. math::
        \mathcal{E} = \frac{1}{2} B \omega l^2
"""


@validate_input(magnetic_induction_=magnetic_field,
    rotation_frequency_=angular_frequency,
    rod_length_=length)
@validate_output(electromotive_force)
def calculate_voltage(magnetic_induction_: Quantity, rotation_frequency_: Quantity,
    rod_length_: Quantity) -> Quantity:
    result_expr = solve(law, electromotive_force, dict=True)[0][electromotive_force]
    result_expr = result_expr.subs({
        magnetic_field: magnetic_induction_,
        angular_frequency: rotation_frequency_,
        length: rod_length_
    })
    return Quantity(result_expr)
