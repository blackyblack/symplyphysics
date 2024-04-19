from sympy import Eq, Derivative, Point2D
from symplyphysics import (
    symbols,
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
## The chemical potential of the system is the amount of energy the system absorbs or releases
## due to the introduction of a particle into the system, i.e. when the particle count increases
## by one.

# Law: mu = (dF/dN)_(T, V)
## mu - chemical potential
## F = F(T, V, N) - [Helmholtz free energy](./helmholtz_free_energy_via_internal_energy.py)
## N - particle count
## T - temperature
## V - volume
## (d/dN)_(T, V) - derivative w.r.t. particle count at constant temperature and volume

chemical_potential = Symbol("chemical_potential", units.energy)
free_energy = Function("free_energy", units.energy)
particle_count = Symbol("particle_count", dimensionless)
temperature = symbols.thermodynamics.temperature
volume = Symbol("volume", units.volume)

# Note that `temperature` and `volume` are only used to show that the free_energy
# function depends on them and they are held constant during differentiation with
# respect to `particle_count`.

law = Eq(
    chemical_potential,
    Derivative(free_energy(temperature, volume, particle_count), particle_count),
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    particle_count_before_=particle_count,
    particle_count_after_=particle_count,
    free_energy_before_=free_energy,
    free_energy_after_=free_energy,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    particle_count_before_: int,
    particle_count_after_: int,
    free_energy_before_: Quantity,
    free_energy_after_: Quantity,
) -> Quantity:
    free_energy_ = two_point_function(
        Point2D(particle_count_before_, free_energy_before_),
        Point2D(particle_count_after_, free_energy_after_),
        particle_count,
    )

    result = law.rhs.subs(
        free_energy(temperature, volume, particle_count),
        free_energy_,
    ).doit()

    return Quantity(result)
