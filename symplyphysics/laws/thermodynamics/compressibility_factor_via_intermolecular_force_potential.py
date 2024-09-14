r"""
Compressibility factor via intermolecular force potential
=========================================================

The virial equation describes the deviation of a real gas from ideal gas behaviour. The virial coefficients
in the virial expansion account for interactions between successively larger groups of molecules. Since
interactions between large numbers of molecules are rare, the virial equation is usually truncated at the third
term onwards. Under the assumption that only pair interactions are present, the compressibility factor can be
linked to the intermolecular force potential.

**Notation:**

#. :math:`k_\text{B}` is the Boltzmann constant.

**Conditions:**

#. The virial expansion is done up to the second virial coefficient inclusively.

..
    TODO Simplify this law by reducing it to the formula of the second virial coefficient
"""

from sympy import Eq, Integral, pi, exp, S
from symplyphysics import (
    clone_as_symbol,
    symbols,
    units,
    dimensionless,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
    convert_to_float,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.chemistry.potential_energy_models import hard_spheres_potential

compressibility_factor = Symbol("compressibility_factor", dimensionless, positive=True)
"""
:doc:`Compressibility factor <definitions.compressibility_factor_is_deviation_from_ideal_gas>` of the gas.

Symbol:
    :code:`Z`
"""

particle_count = Symbol("particle_count", dimensionless, integer=True, positive=True)
"""
Number of particles in the system.

Symbol:
    :code:`N`
"""

volume = Symbol("volume", units.volume, positive=True)
"""
Volume of the system.

Symbol:
    :code:`V`
"""

intermolecular_distance = Symbol("intermolecular_distance", units.length, positive=True)
"""
Distance between gas molecules.

Symbol:
    :code:`r`
"""

intermolecular_force_potential = Function("intermolecular_force_potential", units.energy, real=True)
r"""
Intermolecular force potential as a function of intermolecular distance.

Symbol:
    :code:`phi(r)`

Latex:
    :math:`\varphi(r)`
"""

temperature = clone_as_symbol(symbols.temperature, positive=True)
"""
:attr:`~symplyphysics.symbols.temperature` of the system.
"""

law = Eq(
    compressibility_factor, 1 + 2 * pi * particle_count / volume * Integral(
    (1 - exp(-1 * intermolecular_force_potential(intermolecular_distance) /
    (units.boltzmann_constant * temperature))) * intermolecular_distance**2,
    (intermolecular_distance, 0, S.Infinity),
    ))
r"""
:code:`Z = 1 + 2 * pi * (N / V) * Integral((1 - exp(-1 * phi(r) / (k_B * T))) * r^2, (r, 0, Infinity))`

Latex:
    .. math::
        Z = 1 + \frac{2 \pi N}{V}
                \int \limits_0^{\infty} \left( 1 - \exp{\left( - \frac{\varphi(r)}{k_\text{B} T} \right)} \right) r^2 dr
"""

# Calculate the compressibility factor for the model of hard spheres

_sphere_diameter = Symbol("sphere_diameter", units.length, positive=True)

_hard_spheres_potential = hard_spheres_potential.law.rhs.subs({
    hard_spheres_potential.distance: intermolecular_distance,
    hard_spheres_potential.sphere_diameter: _sphere_diameter,
})

_hard_spheres_compressibility_factor = law.rhs.subs({
    intermolecular_force_potential(intermolecular_distance): _hard_spheres_potential,
}).doit()

# Note that the compressibility factor does not depend on temperature in this model
assert expr_equals(_hard_spheres_compressibility_factor.diff(temperature), 0)


@validate_input(
    number_of_particles_=particle_count,
    volume_=volume,
    sphere_diameter_=intermolecular_distance,
)
@validate_output(compressibility_factor)
def calculate_compressibility_factor(
    number_of_particles_: int,
    volume_: Quantity,
    sphere_diameter_: Quantity,
) -> float:
    # Calculate for the model of hard spheres

    result = _hard_spheres_compressibility_factor.subs({
        particle_count: number_of_particles_,
        volume: volume_,
        _sphere_diameter: sphere_diameter_,
    })
    return convert_to_float(result)
