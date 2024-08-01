r"""
Average speed in Maxwell—Boltzmann statistics
=============================================

The *average*, or mean, *speed* is the expected value of the speed distribution of gas particles.

**Notation:**

#. :math:`k_\text{B}` is the Boltzmann constant.

**Conditions:**

#. The gas is in thermal equilibrium with the environment.
#. The gas particles are distributed according to Maxwell—Boltzmann statistics.
"""

from sympy import Eq, sqrt, pi, S, stats, Interval
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import speed_distribution

average_speed = Symbol("average_speed", units.velocity, positive=True)
r"""
Average molecular speed.

Symbol:
    :code:`avg(v)`

Latex:
    :math:`\langle v \rangle`
"""

equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature,
    "equilibrium_temperature",
    positive=True)
"""
Equilibrium temperature of the gas.

Symbol:
    :code:`T`
"""

molecular_mass = clone_symbol(symbols.basic.mass, "molecular_mass", positive=True)
"""
:attr:`~symplyphysics.symbols.basic.mass` of a gas molecule.

Symbol:
    :code:`m`
"""

law = Eq(
    average_speed,
    sqrt(8 * units.boltzmann_constant * equilibrium_temperature / (pi * molecular_mass)),
)
r"""
:code:`avg(v) = sqrt(8 * k_B * T / (pi * m))`

Latex:
    .. math::
        \langle v \rangle = \sqrt{\frac{8 k_\text{B} T}{\pi m}}
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
