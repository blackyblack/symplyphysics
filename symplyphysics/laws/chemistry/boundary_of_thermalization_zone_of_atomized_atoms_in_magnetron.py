"""
Boundary of thermalization zone of atomized atoms in magnetron
==============================================================

The atoms of the target material evaporate and move towards the substrate inside the
magnetron. There is a distance at which the energy of the traveling atom will decrease to
the energy of the thermal motion of the atoms of the gas-discharge plasma. This distance
is called the boundary of the thermalization zone. The decrease in energy occurs due to
the collision of the traveling atom with the gas atoms in the magnetron. The free path
length is the distance that the sprayed gas travels between two collisions. The traveling
atom moves towards the substrate in the magnetron. At the same time, it collides with gas
atoms. The number of collisions of a traveling atom is the number of collisions, after
which its energy will be equal to the energy of thermal motion in a gas-discharge plasma.

..
    TODO: find link
    TODO: move to `magnetron` folder?
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

boundary_of_thermalization_zone = symbols.length
"""
Boundary of thermalization zone of traveling atoms in magnetron. See :symbols:`length`.
"""

number_of_collisions_of_atom = symbols.nonnegative_number
"""
Number (:symbols:`nonnegative_number`) of collisions of traveling atom with gas atoms.
"""

free_path_length = symbols.mean_free_path
"""
:symbols:`mean_free_path` of traveling atoms.
"""

law = Eq(boundary_of_thermalization_zone, number_of_collisions_of_atom * free_path_length)
"""
:laws:symbol::

:laws:latex::
"""


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
