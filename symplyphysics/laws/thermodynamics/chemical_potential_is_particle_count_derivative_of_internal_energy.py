from sympy import Eq, Derivative, Point2D
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.geometry.line import two_point_function

# Description
## The chemical potential of a system is the amount of energy the system absorbs or releases
## due to the introduction of a particle into the system, i.e. when the particle count increases
## by one.

# Law: mu = (dU/dN)_(S, V)
## mu - chemical potential
## U = U(S, V, N) - internal energy
## N - particle count
## S - entropy
## V - volume
## (d/dN)_(S, V) - derivative w.r.t. particle count at constant entropy and volume

chemical_potential = Symbol("chemical_potential", units.energy)
internal_energy = Function("internal_energy", units.energy)
particle_count = Symbol("particle_count", dimensionless)
entropy = Symbol("entropy", units.energy / units.temperature)
volume = Symbol("volume", units.volume)

# Note that `entropy` and `volume` are only used to show that the internal energy
# function depends on them and they are held constant during differentiation with
# respect to `particle_count`.

law = Eq(
    chemical_potential,
    Derivative(internal_energy(entropy, volume, particle_count), particle_count),
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    particle_count_before_=particle_count,
    particle_count_after_=particle_count,
    internal_energy_before_=internal_energy,
    internal_energy_after_=internal_energy,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    particle_count_before_: int,
    particle_count_after_: int,
    internal_energy_before_: Quantity,
    internal_energy_after_: Quantity,
) -> Quantity:
    internal_energy_ = two_point_function(
        Point2D(particle_count_before_, internal_energy_before_),
        Point2D(particle_count_after_, internal_energy_after_),
        particle_count,
    )

    result = law.rhs.subs(
        internal_energy(entropy, volume, particle_count),
        internal_energy_,
    ).doit()

    return Quantity(result)
