"""
Compressibility factor via intermolecular force potential
=========================================================

The virial equation describes the deviation of a real gas from ideal gas behaviour. The virial coefficients
in the virial expansion account for interactions between successively larger groups of molecules. Since
interactions between large numbers of molecules are rare, the virial equation is usually truncated at the third
term onwards. Under the assumption that only pair interactions are present, the compressibility factor can be
linked to the intermolecular force potential.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The virial expansion is done up to the second virial coefficient inclusively.

**Links:**

#. `Wikipedia, second formula <https://en.wikipedia.org/wiki/Compressibility_factor#Theoretical_models>`__.

..
    TODO Simplify this law by reducing it to the formula of the second virial coefficient
"""

from sympy import Eq, Integral, pi, exp, S
from symplyphysics import (
    clone_as_symbol,
    clone_as_function,
    symbols,
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.chemistry.potential_energy_models import hard_spheres_potential

compressibility_factor = symbols.compressibility_factor
"""
:symbols:`compressibility_factor` of the gas. Also see :doc:`Compressibility factor
<definitions.compressibility_factor_is_deviation_from_ideal_gas>`.
"""

particle_count = clone_as_symbol(
    symbols.particle_count,
    integer=True,
    positive=True,
)
"""
:symbols:`particle_count` of the system.
"""

volume = clone_as_symbol(symbols.volume, positive=True)
"""
:symbols:`volume` of the system.
"""

intermolecular_distance = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="r",
    display_latex="r",
    positive=True,
)
"""
:symbols:`euclidean_distance` between gas molecules.
"""

intermolecular_force_potential = clone_as_function(
    symbols.potential_energy,
    [intermolecular_distance],
    real=True,
)
"""
Intermolecular force potential as a function of :attr:`~intermolecular_distance`. See :symbols:`potential_energy`.
"""

temperature = clone_as_symbol(symbols.temperature, positive=True)
"""
:symbols:`temperature` of the system.
"""

law = Eq(
    compressibility_factor, 1 + 2 * pi * particle_count / volume * Integral(
    (1 - exp(-1 * intermolecular_force_potential(intermolecular_distance) /
    (quantities.boltzmann_constant * temperature))) * intermolecular_distance**2,
    (intermolecular_distance, 0, S.Infinity),
    ))
"""
:laws:symbol::

:laws:latex::
"""

# Calculate the compressibility factor for the model of hard spheres

_sphere_diameter = clone_as_symbol(symbols.diameter, positive=True)

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
