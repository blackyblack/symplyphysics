"""
Number of collisions of particle with gas in magnetron
======================================================

The atoms of the target material evaporate and move towards the substrate inside the
magnetron. The traveling atom moves towards the substrate in the magnetron. At the same
time, it collides with gas atoms. The number of collisions of a traveling atom, after
which its energy will be equal to the energy of thermal motion in a gas-discharge
plasma, can be calculated. This amount will depend on the initial energy of the
traveling atom and the energy transfer coefficient between the atom and the gas atoms.

..
    TODO: find link
    TODO: move to `magnetron` folder?
"""

from sympy import Eq, solve, log
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

collision_count = SymbolNew("N", dimensionless)
"""
Number of collisions between the particle and gas.
"""

initial_energy = clone_as_symbol(symbols.energy, subscript="0")
"""
Initial :symbols:`energy` of the particle.
"""

thermal_energy = symbols.energy
"""
:symbols:`energy` of thermal motion in the plasma.
"""

energy_transfer_coefficient = SymbolNew("x", dimensionless)
"""
Energy transfer coefficient between the particle and the gas.
"""

law = Eq(collision_count,
    log(thermal_energy / initial_energy) / log(1 - energy_transfer_coefficient))
"""
:law:symbol::

:laws:latex::
"""


@validate_input(initial_energy_=initial_energy,
    energy_of_thermal_motion_=thermal_energy,
    energy_transfer_coefficient_=energy_transfer_coefficient)
@validate_output(collision_count)
def calculate_number_of_collisions_of_atoms(initial_energy_: Quantity,
    energy_of_thermal_motion_: Quantity, energy_transfer_coefficient_: float) -> float:
    if energy_transfer_coefficient_ >= 1:
        raise ValueError("Energy transfer coefficient must be less than 1")
    if initial_energy_.scale_factor < energy_of_thermal_motion_.scale_factor:
        raise ValueError(
            "The initial energy of the atom must be greater than or equal to the thermal energy")

    result_expr = solve(law, collision_count,
        dict=True)[0][collision_count]
    result_expr = result_expr.subs({
        initial_energy: initial_energy_,
        thermal_energy: energy_of_thermal_motion_,
        energy_transfer_coefficient: energy_transfer_coefficient_,
    })
    return convert_to_float(result_expr)
