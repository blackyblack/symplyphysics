from sympy import Eq, sqrt, Derivative, I
from sympy.physics.units import hbar
from symplyphysics import (
    units,
    Symbol,
    Function,
    symbols,
    clone_symbol,
)

# Description
## ...

# Notes
## - This law works in the case of 1 spatial dimension. To use it for the 3-dimensional space
##   replace the second derivative with the Laplace operator.

wave_function = Function("wave_function", 1 / sqrt(units.length))
potential_energy = Function("potential_energy", units.energy)
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass")
position = Symbol("position", units.length)
time = Symbol("time", units.time)

law = Eq(
    -1 * hbar**2 / (2 * particle_mass) * Derivative(wave_function(position, time), position, 2)
    + potential_energy(position, time) * wave_function(position, time),
    I * hbar * Derivative(wave_function(position, time), time)
)
