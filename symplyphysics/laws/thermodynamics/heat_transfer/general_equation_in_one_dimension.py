r"""
General equation in one dimension
=================================

Heat equation governs heat diffusion, as well as other diffusive processes. It describes
the evolution of heat transferred from hotter to colder environments in time and space.

**Notes:**

#. There is no straghtforward solution to this equation, and it depends on initial
   conditions as well.
#. To get a similar equation for the 3-dimensional case, replace the spatial
   derivative with gradient :math:`\nabla`.
#. Thermal conductivity :math:`k` can depend not only on position, but also on local
   temperature, but this is out of the scope of this law.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Heat_equation#Non-uniform_isotropic_medium>`__.
"""

from sympy import Eq, Derivative
from symplyphysics import units, Symbol, symbols, clone_as_function

position = symbols.position
"""
:symbols:`position`, or spatial variable.
"""

time = symbols.time
"""
:symbols:`time`.
"""

temperature = clone_as_function(symbols.temperature, [position, time])
"""
:symbols:`temperature` as a function of :attr:`~position` and :attr:`~time`.
"""

medium_density = symbols.density
"""
:symbols:`density` of the medium.
"""

medium_specific_isobaric_heat_capacity = Symbol("c_p", units.energy / (units.temperature * units.mass))
"""
:symbols:`heat_capacity` of the medium at constant :symbols:`pressure` per unit :symbols:`mass`.
"""

thermal_conductivity = clone_as_function(symbols.thermal_conductivity, [position])
"""
:symbols:`thermal_conductivity` of the medium as a function of :attr:`~position`.
"""

heat_source_density = clone_as_function(symbols.energy_density, [position, time], display_symbol="q", display_latex="q")
"""
Density of the rate of heat production by external sources as a function of
:attr:`~position` and :attr:`~time`. See :symbols:`energy_density`.
"""

law = Eq(
    medium_density * medium_specific_isobaric_heat_capacity *
    Derivative(temperature(position, time), time),
    Derivative(
    thermal_conductivity(position) * Derivative(temperature(position, time), position), position) +
    heat_source_density(position, time))
"""
:laws:symbol::

:laws:latex::
"""
