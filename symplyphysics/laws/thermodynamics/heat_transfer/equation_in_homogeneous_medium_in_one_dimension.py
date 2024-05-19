from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Symbol,
    Function,
)

# Description
## Heat equation governs heat diffusion, as well as other diffusive processes. It describes
## the evolution of heat transferred from hotter to colder environments in time and space.

# Law: dT/dt = chi * d^2(T)/dx^2
## T = T(x, t) - temperature
## x - spatial coordinate (position)
## t - time
## chi - [thermal diffusivity](https://en.wikipedia.org/wiki/Thermal_diffusivity)

temperature_function = Function("temperature_function", units.temperature)
position = Symbol("position", units.length)
time = Symbol("time", units.time)
thermal_diffusivity = Symbol("thermal_diffusivity", units.area / units.time)

law = Eq(Derivative(temperature_function(position, time), time),
    thermal_diffusivity * Derivative(temperature_function(position, time), position, 2))

# There is no simple solution to this equation, and it depends on the boundary conditions
# as well.
