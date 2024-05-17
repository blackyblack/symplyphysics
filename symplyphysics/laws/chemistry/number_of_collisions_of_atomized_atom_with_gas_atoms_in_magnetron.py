from sympy import Eq, solve, log
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, convert_to_float)

# Description
## The atoms of the target material evaporate and move towards the substrate inside the magnetron.
## The traveling atom moves towards the substrate in the magnetron. At the same time, it collides with gas atoms.
## The number of collisions of a traveling atom, after which its energy will be equal to the energy of thermal motion
## in a gas-discharge plasma, can be calculated. This amount will depend on the initial energy of the traveling atom and
## the energy transfer coefficient between the atom and the gas atoms.

## Law is: N = ln(Et / E0) / ln(1 - x), where
## N - the number of collisions of traveling atom with gas atoms (this is a statistical quantity, it can be integer, but it can be fractional as well),
## E0 - initial energy of the traveling atom,
## Et - energy of thermal motion in a gas-discharge plasma,
## x - energy transfer coefficient between the traveling atom and the gas atoms.

number_of_collisions_of_atoms = Symbol("number_of_collisions_of_atoms", dimensionless)

initial_energy = Symbol("initial_energy", units.energy)
energy_of_thermal_motion = Symbol("energy_of_thermal_motion", units.energy)
energy_transfer_coefficient = Symbol("energy transfer coefficient", dimensionless)

law = Eq(number_of_collisions_of_atoms,
    log(energy_of_thermal_motion / initial_energy) / log(1 - energy_transfer_coefficient))


def print_law() -> str:
    return print_expression(law)


@validate_input(initial_energy_=initial_energy,
    energy_of_thermal_motion_=energy_of_thermal_motion,
    energy_transfer_coefficient_=energy_transfer_coefficient)
@validate_output(number_of_collisions_of_atoms)
def calculate_number_of_collisions_of_atoms(initial_energy_: Quantity,
    energy_of_thermal_motion_: Quantity, energy_transfer_coefficient_: float) -> float:
    if energy_transfer_coefficient_ >= 1:
        raise ValueError("Energy transfer coefficient must be less than 1")
    if initial_energy_.scale_factor < energy_of_thermal_motion_.scale_factor:
        raise ValueError(
            "The initial energy of the atom must be greater than or equal to the thermal energy")

    result_expr = solve(law, number_of_collisions_of_atoms,
        dict=True)[0][number_of_collisions_of_atoms]
    result_expr = result_expr.subs({
        initial_energy: initial_energy_,
        energy_of_thermal_motion: energy_of_thermal_motion_,
        energy_transfer_coefficient: energy_transfer_coefficient_,
    })
    return convert_to_float(result_expr)
