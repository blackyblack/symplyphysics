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
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
    clone_as_symbol,
)

speed_of_sound = clone_as_symbol(
    symbols.speed,
    display_symbol="v_s",
    display_latex="v_\\text{s}",
)
"""
:symbols:`speed` of sound in the fluid.
"""

density = symbols.density
"""
:symbols:`density` of the fluid.
"""

entropy = symbols.entropy
"""
:symbols:`entropy` of the fluid.
"""

pressure = clone_as_function(symbols.pressure, [density, entropy])
"""
:symbols:`pressure` inside the fluid as a function of :attr:`~density` and
:attr:`~entropy`.
"""

law = Eq(speed_of_sound, sqrt(Derivative(pressure(density, entropy), density)))
"""
:laws:symbol::

:laws:latex::
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
