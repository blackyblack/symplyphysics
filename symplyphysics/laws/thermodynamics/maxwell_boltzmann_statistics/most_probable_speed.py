r"""
Most probable speed in Maxwell—Boltzmann statistics
===================================================

The *most probable speed* is the speed at which the Maxwell—Boltzmann speed distribution function
is maximum.

**Notation:**

#. :math:`k_\text{B}` is the Boltzmann constant.

**Conditions:**

#. The gas is in thermal equilibrium with the environment.
#. The gas particles are distributed according to Maxwell—Boltzmann statistics.
"""

from sympy import Eq, sqrt, solve, sign, S
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

most_probable_speed = Symbol("most_probable_speed", units.velocity, positive=True)
r"""
Most probable speed of particles.

Symbol:
    :code:`v_prob`

Latex:
    :math:`v_\text{prob}`
"""

equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature,
    "equilibrium_temperature",
    positive=True)
"""
Equilibrium :attr:`~symplyphysics.symbols.thermodynamics.temperature` of the gas.

Symbol:
    :code:`T`
"""

molecular_mass = clone_symbol(symbols.basic.mass, "molecular_mass", positive=True)
"""
:attr:`~symplyphysics.symbols.basic.mass` of a gas molecule.

Symbol:
    :code:`m`
"""

law = Eq(most_probable_speed,
    sqrt(2 * units.boltzmann_constant * equilibrium_temperature / molecular_mass))
r"""
:code:`v_prob = sqrt(2 * k_B * T / m)`

Latex:
    .. math::
        v_\text{prob} = \sqrt{\frac{2 k_\text{B} T}{m}}
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
