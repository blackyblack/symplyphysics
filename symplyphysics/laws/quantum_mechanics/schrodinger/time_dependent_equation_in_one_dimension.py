"""
Time dependent Schrödinger equation in one dimension
====================================================

The Schrödinger equation is a linear partial differential equation that governs the wave
function of a quantum-mechanical system. This law describes the general case of a time-dependent
potential and a time-dependent wave function.

**Notation:**

#. :quantity_notation:`hbar`.

**Notes:**

#. This law works in the case of a single spatial dimension. To use it for the :math:`3`-dimensional
   space, replace the spatial second derivative with the Laplace operator.

**Links:**

#. Wikipedia <https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation#Separation_of_variables>
"""

from sympy import Eq, Derivative, I
from symplyphysics import symbols, quantities, clone_as_function

position = symbols.position
"""
:symbols:`position`.
"""

time = symbols.time
"""
:symbols:`time`
"""

wave_function = clone_as_function(symbols.wave_function, [position, time])
"""
Time-dependent :symbols:`wave_function` as a function of :attr:`~position` and :attr:`~time`.
"""

potential_energy = clone_as_function(symbols.potential_energy, [position, time])
"""
Time-independent :symbols:`potential_energy` as a function of :attr:`~position`.
"""

particle_mass = symbols.mass
"""
:symbols:`mass` of the quantum particle.
"""

law = Eq(
    -1 * quantities.hbar**2 / (2 * particle_mass) * Derivative(wave_function(position, time), position, 2) +
    potential_energy(position, time) * wave_function(position, time),
    I * quantities.hbar * Derivative(wave_function(position, time), time))
"""
:laws:symbol::

:laws:latex::
"""
