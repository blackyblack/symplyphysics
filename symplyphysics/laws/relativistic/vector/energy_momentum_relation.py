from sympy import Expr
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    Vector,
    scale_vector,
    QuantityVector,
    vector_magnitude,
    symbols,
)
from symplyphysics.quantities import speed_of_light

# Description
## Relativistic momentum and total relativistic energy are related by a simple equation.

# Law: p * c**2 = E * v
## p - vector of relativistic momentum
## c - speed of light
## E - total energy of relativistic particle
## v - vector of particle's velocity

# Links: Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Energy%E2%80%93momentum_relation#Heuristic_approach_for_massive_particles>
# TODO: find a more exact link

total_energy = symbols.energy


def momentum_law(velocity_: Vector) -> Vector:
    return scale_vector(total_energy / speed_of_light**2, velocity_)


def velocity_law(momentum_: Vector) -> Vector:
    return scale_vector(speed_of_light**2 / total_energy, momentum_)


def total_energy_law(momentum_: Vector, velocity_: Vector) -> Expr:
    p = vector_magnitude(momentum_)
    v = vector_magnitude(velocity_)
    return p * speed_of_light**2 / v


@validate_input(
    total_energy_=total_energy,
    velocity_=units.velocity,
)
@validate_output(units.momentum)
def calculate_momentum(total_energy_: Quantity, velocity_: QuantityVector) -> QuantityVector:
    result_vector = momentum_law(velocity_.to_base_vector())
    return QuantityVector.from_base_vector(
        result_vector,
        subs={total_energy: total_energy_},
    )
