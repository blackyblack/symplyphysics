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

# Law: mu = (dG/dN)_(T, p)
## mu - chemical potential
## G = G(T, p, N) - [Gibbs energy](./isobaric_reaction_potential.py)
## N - particle count
## T - temperature
## p - pressure
## (d/dN)_(T, p) - derivative w.r.t. particle count at constant temperature and pressure

chemical_potential = Symbol("chemical_potential", units.energy)
gibbs_energy = Function("gibbs_energy", units.energy)
particle_count = Symbol("particle_count", dimensionless)
temperature = symbols.thermodynamics.temperature
pressure = Symbol("pressure", units.pressure)

# Note that `temperature` and `pressure` are only used to show that the gibbs_energy
# function depends on them and they are held constant during differentiation with
# respect to `particle_count`.

law = Eq(
    chemical_potential,
    Derivative(gibbs_energy(temperature, pressure, particle_count), particle_count),
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    particle_count_before_=particle_count,
    particle_count_after_=particle_count,
    gibbs_energy_before_=gibbs_energy,
    gibbs_energy_after_=gibbs_energy,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    particle_count_before_: int,
    particle_count_after_: int,
    gibbs_energy_before_: Quantity,
    gibbs_energy_after_: Quantity,
) -> Quantity:
    gibbs_energy_ = two_point_function(
        Point2D(particle_count_before_, gibbs_energy_before_),
        Point2D(particle_count_after_, gibbs_energy_after_),
        particle_count,
    )

    result = law.rhs.subs(
        gibbs_energy(temperature, pressure, particle_count),
        gibbs_energy_,
    ).doit()

    return Quantity(result)
