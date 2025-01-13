"""
Quantum harmonic oscillator equation
====================================

The Schrödinger equation for the quantum simple harmonic oscillator governs the wave function
of the quantum oscillator.

**Notation:**

#. :quantity_notation:`hbar`.

**Links:**

#. `Physics LibreTexts, formula 7.6.4 <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_III_-_Optics_and_Modern_Physics_(OpenStax)/07%3A_Quantum_Mechanics/7.06%3A_The_Quantum_Harmonic_Oscillator>`__.
"""

from sympy import Eq, Derivative
from symplyphysics import symbols, clone_as_function, quantities

position = symbols.position
"""
:symbols:`position` of the particle.
"""

wave_function = clone_as_function(symbols.wave_function, [position])
"""
:symbols:`wave_function` of the oscillating particle.
"""

particle_mass = symbols.mass
"""
:symbols:`mass` of the particle.
"""

particle_energy = symbols.energy
"""
:symbols:`energy` of the particle.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the oscillations.
"""

law = Eq(
    (-1 * quantities.hbar**2 / (2 * particle_mass)) * Derivative(wave_function(position), position, 2) +
    (particle_mass * angular_frequency**2 / 2) * position**2 * wave_function(position),
    particle_energy * wave_function(position),
)
"""
:laws:symbol::

:laws:latex::
"""
