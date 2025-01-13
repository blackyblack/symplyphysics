"""
Time independent solution in one dimension
==========================================

The Schrödinger equation is a linear partial differential equation that governs the wave
function of a quantum-mechanical system.

**Notation:**

#. :quantity_notation:`hbar`.

**Condition:**

#. This law works in the case of a single spatial dimension. To use it for the 3-dimensional space
   replace the spatial second derivative with the Laplace operator.
#. The wave function is independent of time.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation#Separation_of_variables>`__.
"""

from sympy import Eq, Derivative
from symplyphysics import symbols, quantities, clone_as_function

position = symbols.position
"""
:symbols:`position`.
"""

wave_function = clone_as_function(symbols.wave_function, [position])
"""
:symbols:`wave_function` as a function of :attr:`~position`.
"""

potential_energy = clone_as_function(symbols.potential_energy, [position])
"""
:symbols:`potential_energy` as a function of :attr:`~position`.
"""

particle_mass = symbols.mass
"""
:symbols:`mass`.
"""

particle_energy = symbols.energy
"""
:symbols:`energy` of the particle.
"""

law = Eq(
    (-1 * quantities.hbar**2 / (2 * particle_mass)) * Derivative(wave_function(position), position, 2) +
    potential_energy(position) * wave_function(position),
    particle_energy * wave_function(position),
)
"""
:laws:symbol::

:laws:latex::
"""

# The solutions of the Schrödinger equation vary drastically with the kind of potential energy
# function used and also depend on the boundary conditions.
