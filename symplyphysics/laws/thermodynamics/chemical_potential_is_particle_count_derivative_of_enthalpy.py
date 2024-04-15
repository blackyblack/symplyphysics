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
## The chemical potential of the system is the amount of energy the system absorbs or releases
## due to the introduction of a particle into the system, i.e. when the particle count increases
## by one.

# Law: mu = (dH/dN)_(S, p)
## mu - chemical potential
## H = H(S, p, N) - [enthalpy](./enthalpy_is_internal_energy_plus_pressure_energy.py)
## N - particle count
## S - entropy
## p - pressure
## (d/dN)_(S, p) - derivative w.r.t. particle count at constant entropy and pressure

chemical_potential = Symbol("chemical_potential", units.energy)
enthalpy = Function("enthalpy", units.energy)
particle_count = Symbol("particle_count", dimensionless)
entropy = Symbol("entropy", units.energy / units.temperature)
pressure = Symbol("pressure", units.pressure)

# Note that `entropy` and `pressure` are only used to show that the enthalpy
# function depends on them and they are held constant during differentiation with
# respect to `particle_count`.

law = Eq(
    chemical_potential,
    Derivative(enthalpy(entropy, pressure, particle_count), particle_count),
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    particle_count_before_=particle_count,
    particle_count_after_=particle_count,
    enthalpy_before_=enthalpy,
    enthalpy_after_=enthalpy,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    particle_count_before_: int,
    particle_count_after_: int,
    enthalpy_before_: Quantity,
    enthalpy_after_: Quantity,
) -> Quantity:
    enthalpy_ = two_point_function(
        Point2D(particle_count_before_, enthalpy_before_),
        Point2D(particle_count_after_, enthalpy_after_),
        particle_count,
    )

    result = law.rhs.subs(
        enthalpy(entropy, pressure, particle_count),
        enthalpy_,
    ).doit()

    return Quantity(result)
