r"""
Speed distribution
==================

For a system containing a large number of identical non-interacting non-relativistic classical
particles in thermodynamic equilibrium, the speed distribution function is a function such that
:math:`f(v) dv` gives the fraction of particles with speeds in the interval :math:`dv` at speed
:math:`v`.

**Notation:**

#. :math:`k_\text{B}` (:code:`k_B`) is the Boltzmann constant.

**Notes:**

#. Number of particles is big enough that the laws of thermodynamics can be applied.
#. Particles are identical, non-interacting, non-relativistic, and classical.
#. The ensemble of particles is at thermodynamic equilibrium.
"""

from sympy import (
    Eq,
    Rational,
    sqrt,
    pi,
    exp,
    symbols as sym_symbols,
    solve,
    sympify,
)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    clone_symbol,
    symbols,
    CoordinateSystem,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.geometry.elements import volume_element_magnitude
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import velocity_component_distribution

speed_distribution_function = Symbol("speed_distribution_function",
    1 / units.velocity,
    positive=True)
"""
Speed distribution function.

Symbol:
    :code:`f(v)`
"""

particle_speed = Symbol("particle_speed", units.velocity, positive=True)
"""
Particle speed.

Symbol:
    :code:`v`
"""

particle_mass = clone_symbol(symbols.basic.mass, "particle_mass", positive=True)
"""
:attr:`~symplyphysics.symbols.basic.mass` of a particle.

Symbol:
    :code:`m`
"""

equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature,
    "equilibrium_temperature",
    positive=True)
"""
Equilibrium :attr:`~symplyphysics.symbols.thermodynamics.temperature` of the ensemble.

Symbol:
    :code:`T`
"""

law = Eq(
    speed_distribution_function,
    sqrt(2 / pi) * (particle_mass /
    (units.boltzmann_constant * equilibrium_temperature))**Rational(3, 2) * particle_speed**2 *
    exp(-1 * particle_mass * particle_speed**2 /
    (2 * units.boltzmann_constant * equilibrium_temperature)))
r"""
:code:`f(v) = sqrt(2 / pi) * (m / (k_B * T))^(3/2) * v^2 * exp(-1 * m * v^2 / (2 * k_B * T))`

Latex:
    .. math::
        f(v) = \sqrt{\frac{2}{\pi}} \left( \frac{m}{k_\text{B} T} \right)^{3/2} v^2 \exp \left( - \frac{m v^2}{2 k_\text{B} T} \right)
"""

# Derive from Maxwell-Boltzmann distribution of velocity vector components

_velocity_component_distribution = velocity_component_distribution.law.rhs.subs({
    velocity_component_distribution.particle_mass: particle_mass,
    velocity_component_distribution.equilibrium_temperature: equilibrium_temperature,
})

_velocity_x, _velocity_y, _velocity_z = sym_symbols("velocity_x:z", real=True)

_velocity_x_distribution = _velocity_component_distribution.subs(
    velocity_component_distribution.velocity_component, _velocity_x)

_velocity_y_distribution = _velocity_component_distribution.subs(
    velocity_component_distribution.velocity_component, _velocity_y)

_velocity_z_distribution = _velocity_component_distribution.subs(
    velocity_component_distribution.velocity_component, _velocity_z)

_spherical_coordinate_system = CoordinateSystem(CoordinateSystem.System.SPHERICAL)

_radius, _azimuthal_angle, _polar_angle = _spherical_coordinate_system.coord_system.base_scalars()

# The speed distribution depends solely on the radial component of the velocity vector.
# Therefore we integrate over polar and azimuthal angles to get rid of them.
# We work in the three-dimensional velocity space `d^3(v)` so radius is speed (magnitude of velocity vector)
_spherical_volume_velocity_element = (sympify(
    volume_element_magnitude(_spherical_coordinate_system)).subs(_radius, particle_speed).integrate(
    (_polar_angle, 0, pi), (_azimuthal_angle, 0, 2 * pi)))

# `v_x`, `v_y` and `v_z` are independent random variables, therefore we can get the distribution of
# the velocity vector `v` by multiplying the distributions of its coordinates.
_cartesian_speed_distribution = _velocity_x_distribution * _velocity_y_distribution * _velocity_z_distribution

# In order to switch to spherical coordinates we multiply it by the spherical volume element of the
# 3D velocity space found above.
_spherical_speed_distribution = _cartesian_speed_distribution * _spherical_volume_velocity_element

_speed_eqn = Eq(particle_speed**2, _velocity_x**2 + _velocity_y**2 + _velocity_z**2)

_speed_distribution_derived = solve(
    [
    Eq(speed_distribution_function, _spherical_speed_distribution),
    _speed_eqn,
    ],
    (speed_distribution_function, _velocity_x),
    dict=True,
)[0][speed_distribution_function]

assert expr_equals(_speed_distribution_derived, law.rhs)


@validate_input(
    particle_speed_=particle_speed,
    particle_mass_=particle_mass,
    ensemble_temperature_=equilibrium_temperature,
)
@validate_output(speed_distribution_function)
def calculate_speed_distribution_function(
    particle_speed_: Quantity,
    particle_mass_: Quantity,
    ensemble_temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        particle_speed: particle_speed_,
        particle_mass: particle_mass_,
        equilibrium_temperature: ensemble_temperature_,
    })
    return Quantity(result)
