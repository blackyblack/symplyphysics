"""
Magnetic field due to finite coil along axis
============================================

Using the Biot—Savart law, it is possible to obtain the formula for the magnetic flux density at
any point on the axis of a coil. It is directly proportional to the coil's turn count and current
and inversely proportional to its length, and also depends on the position of the measurement
point relative to the coil ends.

**Conditions:**

#.  The point of measurement must lie on the axis of the coil.

#.  The medium is a vacuum.

#.  The :math:`z`-axis is the axis of rotation and is oriented according to the right-hand side
    rule.

**Links:**

#. `Physics LibreTexts — Solenoids and Toroids. Equation (12.7.4) <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/12%3A_Sources_of_Magnetic_Fields/12.07%3A_Solenoids_and_Toroids>`__.
"""

from sympy import (Eq, cos)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, quantities,
    clone_as_symbol)
from symplyphysics.core.symbols.quantities import scale_factor

magnetic_flux_density = symbols.magnetic_flux_density
"""
Magnitude of :symbols:`magnetic_flux_density`.
"""

current = symbols.current
"""
:symbols:`current` flowing through the wire.
"""

turn_count = symbols.positive_number
"""
Number of turns in the coil. See :symbols:`positive_number`.
"""

coil_length = clone_as_symbol(symbols.length, display_latex="\\ell")
"""
:symbols:`length` of the coil.
"""

first_angle = clone_as_symbol(symbols.angle, subscript="1")
"""
Acute :symbols:`angle` between the coil axis (or side) and the vector from the measuring point and
the first end of the coil (that has a smaller :math:`z` coordinate).
"""

second_angle = clone_as_symbol(symbols.angle, subscript="2")
"""
Acute :symbols:`angle` between the coil axis (or side) and the vector from the measuring point and
the second end of the coil (that has a greater :math:`z` coordinate).
"""

law = Eq(
    magnetic_flux_density,
    (quantities.vacuum_permeability * current * turn_count / (2 * coil_length)) *
    (cos(first_angle) + cos(second_angle)),
)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive from the magnetic field due to a current loop.


@validate_input(
    current_=current,
    turn_count_=turn_count,
    coil_length_=coil_length,
    first_angle_=first_angle,
    second_angle_=second_angle, 
)
@validate_output(magnetic_flux_density)
def calculate_magnetic_flux_density(
    current_: Quantity,
    turn_count_: int,
    coil_length_: Quantity,
    first_angle_: float | Quantity,
    second_angle_: float | Quantity,
) -> Quantity:
    result = law.rhs.subs({
        current: current_,
        turn_count: turn_count_,
        coil_length: coil_length_,
        first_angle: scale_factor(first_angle_),
        second_angle: scale_factor(second_angle_),
    })

    return Quantity(result)
