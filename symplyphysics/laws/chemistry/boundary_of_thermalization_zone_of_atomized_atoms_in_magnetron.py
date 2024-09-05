from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
)

# Description
## The atoms of the target material evaporate and move towards the substrate inside the magnetron. There is a distance at which the energy of the traveling
## atom will decrease to the energy of the thermal motion of the atoms of the gas-discharge plasma. This distance is called the boundary of the thermalization zone.
## The decrease in energy occurs due to the collision of the traveling atom with the gas atoms in the magnetron. The free path length is the distance that the sprayed
## gas travels between two collisions.
## The traveling atom moves towards the substrate in the magnetron. At the same time, it collides with gas atoms.
## The number of collisions of a traveling atom is the number of collisions, after which its energy will be equal to the energy of thermal motion
## in a gas-discharge plasma.

## Law is: L = N * l, where
## L - boundary of thermalization zone of traveling atoms in magnetron,
## N - the number of collisions of traveling atom with gas atoms (this is a statistical quantity, it can be integer, but it can be fractional as well),
## l - free path length of traveling atom.

boundary_of_thermalization_zone = Symbol("boundary_of_thermalization_zone", units.length)

number_of_collisions_of_atom = Symbol("number_of_collisions_of_atom", dimensionless)
free_path_length = Symbol("free_path_length", units.length)

law = Eq(boundary_of_thermalization_zone, number_of_collisions_of_atom * free_path_length)


@validate_input(number_of_collisions_of_atom_=number_of_collisions_of_atom,
    free_path_length_=free_path_length)
@validate_output(boundary_of_thermalization_zone)
def calculate_boundary_of_thermalization_zone(number_of_collisions_of_atom_: float,
    free_path_length_: Quantity) -> Quantity:
    result_expr = solve(law, boundary_of_thermalization_zone,
        dict=True)[0][boundary_of_thermalization_zone]
    result_expr = result_expr.subs({
        number_of_collisions_of_atom: number_of_collisions_of_atom_,
        free_path_length: free_path_length_,
    })
    return Quantity(result_expr)
