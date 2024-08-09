"""
Isentropic speed of sound
=========================

Derived by Laplace, the formula for the speed of sound in fluids uses the pressure-density dependence
and the fact that the oscillations in a sound wave happen so fast, and the thermal conductivity of air
is so small, that there is no heat transfer in the sound wave, i.e. it is an adiabatic, and therefore
isentropic, process. This is in contrast with the Newton's formula, who thought that the sound propagation
is an isothermal process in the assumption that the temperature differences between different parts of
the sound wave immediately level out, which eventually turned out to be inconsistent with experimental data.
"""

from sympy import Eq, sqrt, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)

speed_of_sound = Symbol("speed_of_sound", units.velocity)
"""
Speed of sound in the fluid.

Symbol:
    :code:`c`
"""

pressure = Function("pressure", units.pressure)
"""
Pressure inside the fluid.

Symbol:
    :code:`p`
"""

density = Symbol("density", units.mass / units.volume)
r"""
Density of the fluid.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

entropy = Symbol("entropy", units.energy / units.temperature)
"""
Entropy of the fluid.

Symbol:
    :code:`S`
"""

law = Eq(speed_of_sound, sqrt(Derivative(pressure(density, entropy), density)))
r"""
:code:`c = sqrt(Derivative(p(rho, S), rho))`

Latex:
    .. math::
        c = \sqrt{\left( \frac{\partial p}{\partial \rho} \right)_S}
"""


@validate_input(
    density_change_=density,
    pressure_change_=pressure,
)
@validate_output(speed_of_sound)
def calculate_speed_of_sound(
    density_change_: Quantity,
    pressure_change_: Quantity,
) -> Quantity:
    pressure_ = (pressure_change_ / density_change_) * density
    result = law.rhs.subs(pressure(density, entropy), pressure_).doit()
    return Quantity(result)
