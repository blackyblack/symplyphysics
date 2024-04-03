from sympy import Eq, Integral, pi, exp, S, Piecewise
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    convert_to,
)
from symplyphysics.core.expr_comparisons import expr_equals

# Description
## The virial equation describes the deviation of a real gas from ideal gas behaviour. The virial coefficients
## in the virial expansion account for interactions between successively larger groups of molecules. Since
## interactions between large numbers of molecules are rare, the virial equation is usually truncated at the third
## term onwards. Under the assumption that only pair interactions are present, the compressibility factor can be
## linked to the intermolecular force potential:

# Law: Z = 1 + 2 * pi * N / V * Integral((1 - exp(-phi(r) / (k * T))) * r**2, (r, 0, infinity))
## Z - [compressibility factor](../../definitions/compressibility_factor_is_deviation_from_ideal_gas.py)
## N - number of particles
## V - volume
## r - radius
## phi(r) - intermolecular force potential as a function of radius
## k - Boltzmann constant
## T - absolute temperature

# Conditions
## - The virial expansion is done up to the second virial coefficient inclusively.

compressibility_factor = Symbol("compressibility_factor", dimensionless, positive=True)
number_of_particles = Symbol("number_of_particles", dimensionless, integer=True, positive=True)
volume = Symbol("volume", units.volume, positive=True)
radius = Symbol("radius", units.length, positive=True)
intermolecular_force_potential = Function("intermolecular_force_potential", units.energy, real=True)
temperature = Symbol("temperature", units.temperature, positive=True)

law = Eq(
    compressibility_factor,
    1 + 2 * pi * number_of_particles / volume
    * Integral(
        (1 - exp(-1 * intermolecular_force_potential(radius) / (units.boltzmann_constant * temperature))) * radius**2,
        (radius, 0, S.Infinity),
    )
)

# Calculate the compressibility factor for the model of hard spheres

_sphere_radius = Symbol("sphere_radius", units.length, positive=True)

_hard_spheres_potential = Piecewise(
    (S.Infinity, radius < _sphere_radius),
    (0, radius >= _sphere_radius),
)

_hard_spheres_compressibility_factor = law.rhs.subs(
    intermolecular_force_potential(radius), _hard_spheres_potential
).doit()

# Notice that the compressibility factor does not depend on temperature in this model
assert expr_equals(_hard_spheres_compressibility_factor.diff(temperature), 0)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    number_of_particles_=number_of_particles,
    volume_=volume,
    sphere_radius_=radius,
)
@validate_output(compressibility_factor)
def calculate_compressibility_factor(
    number_of_particles_: int | float,  # float to allow numbers written in exponent notation
    volume_: Quantity,
    sphere_radius_: Quantity,
) -> float:
    # Calculate for the model of hard spheres

    result = _hard_spheres_compressibility_factor.subs({
        number_of_particles: number_of_particles_,
        volume: volume_,
        _sphere_radius: sphere_radius_,
    })
    return float(convert_to(Quantity(result), S.One))
