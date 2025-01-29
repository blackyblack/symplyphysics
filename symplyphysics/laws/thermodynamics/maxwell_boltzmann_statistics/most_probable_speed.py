r"""
Most probable speed in Maxwell—Boltzmann statistics
===================================================

The *most probable speed* is the speed at which the Maxwell—Boltzmann speed distribution function
is maximum.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The gas is in thermal equilibrium with the environment.
#. The gas particles are distributed according to Maxwell—Boltzmann statistics.

**Links:**

#. `Wikipedia, first item within the list <https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Typical_speeds>`__.
"""

from sympy import Eq, sqrt, solve, sign, S
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import speed_distribution

most_probable_speed = clone_as_symbol(symbols.speed, subscript="\\text{prob}")
"""
Most probable :symbols:`speed` of particles.
"""

equilibrium_temperature = clone_as_symbol(symbols.temperature, positive=True)
"""
Equilibrium :symbols:`temperature` of the gas.
"""

molecular_mass = clone_as_symbol(symbols.mass, positive=True)
"""
:symbols:`mass` of a gas molecule.
"""

law = Eq(most_probable_speed,
    sqrt(2 * quantities.boltzmann_constant * equilibrium_temperature / molecular_mass))
"""
:laws:symbol::

:laws:latex::
"""

# Derive from the Maxwell-Boltzmann speed distribution function

_distribution = speed_distribution.law.rhs.subs({
    speed_distribution.particle_mass: molecular_mass,
    speed_distribution.equilibrium_temperature: equilibrium_temperature,
})

_distribution_first_derivative = _distribution.diff(speed_distribution.particle_speed)

# Found the points of extremum of the speed distribution function
_solutions = solve(_distribution_first_derivative, speed_distribution.particle_speed)

# Asserting there is only one point of extremum
assert len(_solutions) == 1
_most_probable_speed_derived = _solutions[0]

assert expr_equals(_most_probable_speed_derived, law.rhs)

_distribution_second_derivative = _distribution_first_derivative.diff(
    speed_distribution.particle_speed)

_distribution_second_derivative_at_most_probable_speed = _distribution_second_derivative.subs(
    speed_distribution.particle_speed, _most_probable_speed_derived)

# Proved that the point of extremum found is the point of maximum
assert sign(_distribution_second_derivative_at_most_probable_speed) == S.NegativeOne


@validate_input(
    equilibrium_temperature_=equilibrium_temperature,
    particle_mass_=molecular_mass,
)
@validate_output(most_probable_speed)
def calculate_most_probable_speed(
    equilibrium_temperature_: Quantity,
    particle_mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        equilibrium_temperature: equilibrium_temperature_,
        molecular_mass: particle_mass_,
    })
    return Quantity(result)
