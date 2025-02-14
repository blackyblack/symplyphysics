r"""
Average speed in Maxwell—Boltzmann statistics
=============================================

The *average*, or mean, *speed* is the expected value of the speed distribution of gas particles.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The gas is in thermal equilibrium with the environment.
#. The gas particles are distributed according to Maxwell—Boltzmann statistics.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Typical_speeds>`__.
"""

from sympy import Eq, sqrt, pi, S, stats, Interval
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

average_speed = clone_as_symbol(symbols.speed, display_symbol="avg(v)", display_latex="\\langle v \\rangle")
"""
Average molecular :symbols:`speed`.
"""

equilibrium_temperature = clone_as_symbol(symbols.temperature, positive=True)
"""
Equilibrium :symbols:`temperature` of the gas.
"""

molecular_mass = clone_as_symbol(symbols.mass, positive=True)
"""
:symbols:`mass` of a gas molecule.
"""

law = Eq(
    average_speed,
    sqrt(8 * quantities.boltzmann_constant * equilibrium_temperature / (pi * molecular_mass)),
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from Maxwell-Boltzmann distribution function

_distribution = speed_distribution.law.rhs.subs({
    speed_distribution.equilibrium_temperature: equilibrium_temperature,
    speed_distribution.particle_mass: molecular_mass,
})

_speed_random_variable = stats.ContinuousRV(
    speed_distribution.particle_speed,
    _distribution,
    set=Interval(0, S.Infinity),
)

# Average speed is the expected value of the speed random variable
_average_speed_derived = stats.E(_speed_random_variable)

assert expr_equals(_average_speed_derived, law.rhs)


@validate_input(
    equilibrium_temperature_=equilibrium_temperature,
    particle_mass_=molecular_mass,
)
@validate_output(average_speed)
def calculate_average_speed(
    equilibrium_temperature_: Quantity,
    particle_mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        equilibrium_temperature: equilibrium_temperature_,
        molecular_mass: particle_mass_,
    })
    return Quantity(result)
