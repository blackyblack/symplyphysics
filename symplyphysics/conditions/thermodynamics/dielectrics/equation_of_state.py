"""
Equation of state
=================

In the physics of dielectric materials, the "equation of state" is the equation that connects
the electric displacement with electric field strength, temperature, and density of the
medium. This equation cannot be determined solely from thermodynamics and must be either
supplied from experiment or derived from the theory of dielectric polarization.

**Links:**

#. Formula 31.10 on p. 122 of "General Course of Physics" (Obschiy kurs fiziki), vol. 3 by Sivukhin D.V. (1979).
"""

from sympy import Eq
from symplyphysics import FunctionNew, symbols

electric_displacement = symbols.electric_displacement
"""
:symbols:`electric_displacement`.
"""

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the medium.
"""

density = symbols.density
"""
:symbols:`density` of the medium.
"""

state_function = FunctionNew(
    "f",
    [electric_field_strength, temperature, density],
    electric_displacement.dimension,
)
"""
The function that describes the relationship between :attr:`~electric_displacement` and
:attr:`~electric_field_strength`, :attr:`~temperature` and :attr:`~density`.
"""

condition = Eq(electric_displacement, state_function(electric_field_strength, temperature, density))
"""
:laws:symbol::

:laws:latex::
"""
