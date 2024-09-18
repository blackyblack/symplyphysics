r"""
Electromotive force induced in rotating rod
===========================================

Let a rod rotate in a uniform magnetic field. The plane of rotation is perpendicular
to the magnetic field lines. The axis of rotation passes through one of the ends of
the rod. Then the electromotive force induced at the ends of the rod depends on the
magnitude of the magnetic flux density, the rotation frequency and the length of the rod.

**Conditions:**

#. The angular velocity of the rod is parallel to the magnetic field. This means that
   the rod is rotating in a plane perpendicular to the magnetic field.
#. The magnetic field is uniform.
#. The angular velocity of the rod is constant.

**Links:**

#. `Example 13.4.2 <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/13%3A_Electromagnetic_Induction/13.04%3A_Motional_Emf>`_.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

electromotive_force = symbols.electromotive_force
"""
:symbols:`electromotive_force` induced in the rod.
"""

magnetic_flux_density = symbols.magnetic_flux_density
"""
Magnitude of :symbols:`magnetic_flux_density`.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of rod's rotation.
"""

length = symbols.length
"""
:symbols:`length` of the rod.
"""

law = Eq(electromotive_force, magnetic_flux_density * angular_frequency * length**2 / 2)
"""
:laws:symbol::

:laws:latex::
"""


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
