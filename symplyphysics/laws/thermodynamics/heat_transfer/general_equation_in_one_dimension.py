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
"""

from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Symbol,
    Function,
)

temperature = Function("temperature", units.temperature)
"""
Temperature as a function of position and time.

Symbol:
    :code:`T(x, t)`
"""

medium_density = Symbol("medium_density", units.mass / units.volume)
r"""
Density of the medium.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

medium_specific_isobaric_heat_capacity = Symbol("medium_specific_isobaric_heat_capacity",
    units.energy / (units.temperature * units.mass))
"""
Heat capacity of the medium at constant pressure per unit mass.

Symbol:
    :code:`c_p`

Latex:
    :math:`c_p`
"""

thermal_conductivity = Function(
    "thermal_conductivity",
    units.power / (units.length * units.temperature),
)
"""
`Thermal conductivity <https://en.wikipedia.org/wiki/Thermal_conductivity_and_resistivity#Definition>`_
of the medium as a function of position.

Symbol:
    :code:`k(x)`
"""

heat_source_density = Function(
    "heat_source_density",
    units.power / units.volume,
)
"""
Density of the rate of heat production by external sources as a function of position and time.

Symbol:
    :code:`q(x, t)`
"""

position = Symbol("position", units.length)
"""
Position, or spatial variable.

Symbol:
    :code:`x`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

law = Eq(
    medium_density * medium_specific_isobaric_heat_capacity *
    Derivative(temperature(position, time), time),
    Derivative(
    thermal_conductivity(position) * Derivative(temperature(position, time), position),
    position) + heat_source_density(position, time))
r"""
:code:`rho * c_p * Derivative(T(x, t), t) = Derivative(k(x) * Derivative(T(x, t), x), x) + q(x, t)`

Latex:
    .. math::
        \rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x} \left( k(x) \frac{\partial T}{\partial x} \right) + q(x, t)
"""
