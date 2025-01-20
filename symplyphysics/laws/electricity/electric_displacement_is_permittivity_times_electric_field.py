"""
Electric displacement is permittivity times electric field
==========================================================

In a linear, homogeneous, isotropic dielectric the electric displacement is linearly
proportional to the electric field strength, with the constant of proprortionality
being the permittivity of the medium.

**Notes:**

#. If the medium is anisotropic, the relation between electric displacement and
   electric field is similar, but permittivity is now a tensor and not a scalar.

#. In a nonhomogeneous medium, permittivity is a function of position inside the medium.

#. In a nonlinear medium, permittivity is a function of the electric field and has a
   time-dependent response.

**Conditions:**

#. The medium is linear, homogeneous, and isotropic.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Electric_displacement_field#Definition>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

electric_displacement = symbols.electric_displacement
"""
:symbols:`electric_displacement` in the medium.
"""

absolute_permittivity = symbols.absolute_permittivity
"""
:symbols:`absolute_permittivity` of the medium.
"""

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength` in the medium.
"""

law = Eq(electric_displacement, absolute_permittivity * electric_field_strength)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    absolute_permittivity_=absolute_permittivity,
    electric_field_strength_=electric_field_strength,
)
@validate_output(electric_displacement)
def calculate_electric_displacement(
    absolute_permittivity_: Quantity,
    electric_field_strength_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        absolute_permittivity: absolute_permittivity_,
        electric_field_strength: electric_field_strength_,
    })
    return Quantity(result)
