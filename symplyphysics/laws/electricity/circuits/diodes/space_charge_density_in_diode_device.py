"""
Charge density in diode
=======================

The simplest device implementing a cathode current control method is a diode device. It
contains a vacuum chamber, a thermionic cathode, a mesh electrode and an anode.

**Notation:**

#. :quantity_notation:`vacuum_permittivity`.

..
    TODO: find link
"""

from sympy import Eq, solve, Rational
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    quantities,
    symbols,
)

charge_density = symbols.volumetric_charge_density
"""
:symbols:`volumetric_charge_density` in the diode.
"""

grid_voltage = symbols.voltage
"""
:symbols:`voltage` on the grid.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between the grid and the cathode.
"""

law = Eq(charge_density,
    Rational(4, 9) * quantities.vacuum_permittivity * grid_voltage / distance**2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(voltage_on_grid_=grid_voltage,
    distance_between_cathode_and_grid_=distance)
@validate_output(charge_density)
def calculate_space_charge_density(voltage_on_grid_: Quantity,
    distance_between_cathode_and_grid_: Quantity) -> Quantity:
    result_expr = solve(law, charge_density, dict=True)[0][charge_density]
    result_expr = result_expr.subs({
        grid_voltage: voltage_on_grid_,
        distance: distance_between_cathode_and_grid_,
    })
    return Quantity(result_expr)
