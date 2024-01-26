from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The work-energy principle states that the work done by all forces acting on a particle
## (the work of the resultant force) equals the change in the kinetic energy of the particle.

# Law: W = K_after - K_before
## W - total work on particle
## K_after, K_before - kinetic energy of particle before and after work is done, respectively

total_work = Symbol("total_work", units.energy)
kinetic_energy = Function("kinetic_energy", units.energy)
time_before = Symbol("time_before", units.time)
time_after = Symbol("time_after", units.time)

law = Eq(total_work, kinetic_energy(time_after) - kinetic_energy(time_before))


def print_law() -> str:
    return print_expression(law)


@validate_input(kinetic_energy_before_=kinetic_energy, kinetic_energy_after_=kinetic_energy)
@validate_output(total_work)
def calculate_total_work(kinetic_energy_before_: Quantity, kinetic_energy_after_: Quantity) -> Quantity:
    result = law.rhs.subs({
        kinetic_energy(time_before): kinetic_energy_before_,
        kinetic_energy(time_after): kinetic_energy_after_,
    })
    return Quantity(result)
