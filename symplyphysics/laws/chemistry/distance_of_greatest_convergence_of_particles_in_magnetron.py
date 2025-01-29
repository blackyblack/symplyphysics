"""
Distance of greatest convergence of particles in magnetron
==========================================================

The traveling atom moves towards the substrate in the magnetron. At the same time, it
collides with gas atoms. The energy transfer coefficient in these collisions depends on
the mass of the traveling atom and the mass of the gas atom. The distance of the
greatest convergence of colliding particles can be calculated using the model of
quasi-rigid spheres. The discharge voltage is the voltage between the cathode and the
anode in the magnetron at which plasma occurs.

..
    TODO: find link
    TODO: move to `magnetron` folder?
"""

from sympy import Eq, solve, log
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

distance_of_convergence = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` of greatest convergence of two colliding particles.
"""

discharge_voltage = symbols.voltage
"""
Discharge :symbols:`voltage`.
"""

first_atomic_number = clone_as_symbol(symbols.atomic_number, subscript="1")
"""
:symbols:`atomic_number` of first atom.
"""

second_atomic_number = clone_as_symbol(symbols.atomic_number, subscript="2")
"""
:symbols:`atomic_number` of second atom.
"""

distance_constant = Quantity(0.122e-10 * units.meter, display_symbol="d_0")
"""
Constant equal to :math:`0.122 \\cdot 10^{-10} \\, \\text{m}`.
"""

voltage_constant = Quantity(95.863 * units.volt, display_symbol="V_0")
"""
Constant equal to :math:`95.863 \\, \\text{V}`.
"""

law = Eq(
    distance_of_convergence, -1 * distance_constant *
    (first_atomic_number**0.0387 + second_atomic_number**0.0387) * log(discharge_voltage /
    (voltage_constant * (first_atomic_number * second_atomic_number)**1.4883)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(discharge_voltage_=discharge_voltage,
    atomic_number_of_first_atom_=first_atomic_number,
    second_atomic_number_=second_atomic_number)
@validate_output(distance_of_convergence)
def calculate_distance_of_convergence_of_particles(discharge_voltage_: Quantity,
    atomic_number_of_first_atom_: int, second_atomic_number_: int) -> Quantity:
    result_expr = solve(law, distance_of_convergence, dict=True)[0][distance_of_convergence]
    result_expr = result_expr.subs({
        discharge_voltage: discharge_voltage_,
        first_atomic_number: atomic_number_of_first_atom_,
        second_atomic_number: second_atomic_number_,
    })
    return Quantity(result_expr)
