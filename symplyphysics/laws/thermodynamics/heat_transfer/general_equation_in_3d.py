"""
General heat equation in 3D
===========================

Heat equation governs heat diffusion, as well as other diffusive processes. It describes the
evolution of heat transferred from hotter to colder environments in time and space.

**Conditions:**

#. The medium is isotropic.

#. There is no mass transfer or radiation in the system.

**Links:**

#. `Wikipedia — Heat equation <https://en.wikipedia.org/wiki/Heat_equation#Non-uniform_isotropic_medium>`__.
"""

from sympy import Eq, Derivative, evaluate
from symplyphysics import symbols, clone_as_function, Symbol, Function, units

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.operators import VectorDivergence, VectorGradient

medium_density = symbols.density
"""
:symbols:`density` of the medium.
"""

medium_specific_isobaric_heat_capacity = Symbol(
    "c_p",
    units.energy / (units.mass * units.temperature),
)
"""
:symbols:`mass`-specific isobaric :symbols:`heat_capacity` of the medium.
"""

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Position vector. See :symbols:`distance_to_origin`.
"""

time = symbols.time
"""
:symbols:`time`.
"""

temperature = clone_as_function(symbols.temperature, [position_vector, time])
"""
:symbols:`temperature` as a function of :attr:`~position_vector` and :attr:`~time`.
"""

thermal_conductivity = clone_as_function(symbols.thermal_conductivity, [position_vector])
"""
:symbols:`thermal_conductivity` as a function of :attr:`~position_vector`.
"""

heat_source_density = Function("q", [position_vector], dimension=units.power / units.volume)
"""
Volumetric density of heat sources (i.e. :symbols:`power` produced per unit :symbols:`volume`) as
a function of :attr:`~position_vector`.
"""

with evaluate(False):
    _rho = medium_density
    _c_p = medium_specific_isobaric_heat_capacity
    _dt_dt = Derivative(temperature(position_vector, time), time)
    _k = thermal_conductivity(position_vector)
    _grad_t = VectorGradient(temperature(position_vector, time), evaluate=False)

law = Eq(
    _rho * _c_p * _dt_dt,
    VectorDivergence(_k * _grad_t, evaluate=False) + heat_source_density(position_vector),
)
"""
:laws:symbol::

:laws:latex::
"""


# UNIQUE_LAW_ID: 130
