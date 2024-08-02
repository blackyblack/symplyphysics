from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Symbol,
    Function,
)

# Description
## Heat equation governs heat diffusion, as well as other diffusive processes. It describes
## the evolution of heat transferred from hotter to colder environments in time and space.

# Law: rho * c_v * dT/dt = d(k(x) * dT/dx)/dx + q(x, t)
## T = T(x, t) - temperature
## rho - density of medium
## c_v - mass-specific heat capacity of medium
## k - [thermal conductivity](https://en.wikipedia.org/wiki/Thermal_conductivity_and_resistivity#Definition) of medium
## q - density of heat sources
## x - position, spatial variable
## t - time
## d/dt - time derivative
## d/dx - spatial derivative

# Notes
## - To get a similar equation for the 3-dimensional case, replace the spatial
##   derivative with the gradient operator.
## - Thermal conductivity can depend not only on position, but also on local temperature, but this case is out
##   of the scope of this law.

temperature_function = Function("temperature_function", units.temperature)
medium_density = Symbol("medium_density", units.mass / units.volume)
medium_specific_heat_capacity = Symbol("medium_specific_heat_capacity",
    units.energy / (units.temperature * units.mass))
thermal_conductivity = Function(
    "thermal_conductivity",
    units.power / (units.length * units.temperature),
)
heat_source_density = Function(
    "heat_source_density",
    units.power / units.volume,
)
position = Symbol("position", units.length)
time = Symbol("time", units.time)

law = Eq(
    medium_density * medium_specific_heat_capacity *
    Derivative(temperature_function(position, time), time),
    Derivative(
    thermal_conductivity(position) * Derivative(temperature_function(position, time), position),
    position) + heat_source_density(position, time))
