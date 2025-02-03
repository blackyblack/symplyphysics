"""
Equivalent diode voltage for triode
===================================

A triode has three electrodes: a cathode, an anode and one control grid. The triode can
be replaced with an equivalent diode and the :math:`3/2`-power law can be applied. The
voltage in the equivalent diode can also be calculated.

..
    TODO: find link
"""

from sympy import Eq, Rational, solve
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_symbol)

equivalent_diode_voltage = symbols.voltage
"""
:symbols:`voltage` between the cathode and anode of an equivalent diode for a triode.
"""

anode_distance = clone_as_symbol(symbols.euclidean_distance, display_symbol="d_a", display_latex="d_\\text{a}")
"""
:symbols:`euclidean_distance` between the cathode and the anode.
"""

grid_distance = clone_as_symbol(symbols.euclidean_distance, display_symbol="d_g", display_latex="d_\\text{g}")
"""
:symbols:`euclidean_distance` between the cathode and the grid.
"""

anode_voltage = clone_as_symbol(symbols.voltage, display_symbol="U_a", display_latex="U_\\text{a}")
"""
:symbols:`voltage` between the cathode and the anode.
"""

voltage_gain = clone_as_symbol(symbols.circuit_gain, subscript="V")
"""
Voltage gain (:symbols:`circuit_gain`) of the triode.
"""

grid_voltage = clone_as_symbol(symbols.voltage, display_symbol="U_g", display_latex="U_\\text{g}")
"""
:symbols:`voltage` between the cathode and the grid.
"""

law = Eq(equivalent_diode_voltage, (grid_voltage + anode_voltage / voltage_gain) / (1 +
    ((anode_distance / grid_distance)**Rational(4, 3)) / voltage_gain))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(distance_to_anode_=anode_distance,
    distance_to_grid_=grid_distance,
    anode_voltage_=anode_voltage,
    voltage_triode_gain_=voltage_gain,
    grid_voltage_=grid_voltage)
@validate_output(equivalent_diode_voltage)
def calculate_voltage_of_equivalent_diode(distance_to_anode_: Quantity, distance_to_grid_: Quantity,
    anode_voltage_: Quantity, voltage_triode_gain_: float, grid_voltage_: Quantity) -> Quantity:
    result_expr = solve(law, equivalent_diode_voltage, dict=True)[0][equivalent_diode_voltage]
    result_expr = result_expr.subs({
        anode_distance: distance_to_anode_,
        grid_distance: distance_to_grid_,
        anode_voltage: anode_voltage_,
        voltage_gain: voltage_triode_gain_,
        grid_voltage: grid_voltage_,
    })
    return Quantity(result_expr)
