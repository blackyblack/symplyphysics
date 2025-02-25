"""
Voltage of equivalent diode for pentode
=======================================

A pentode has five electrodes: a cathode, an anode, and three grids. The pentode can be
replaced with an equivalent diode, and the voltage in the equivalent diode can be
calculated.

..
    TODO: find link
"""

from sympy import Eq, Rational, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

equivalent_diode_voltage = symbols.voltage
"""
:symbols:`voltage` between the cathode and anode of an equivalent diode for a pentode.
"""

first_grid_voltage = clone_as_symbol(symbols.voltage, subscript="1")
"""
:symbols:`voltage` between the cathode and the first grid.
"""

second_grid_voltage = clone_as_symbol(symbols.voltage, subscript="2")
"""
:symbols:`voltage` between the cathode and the second grid.
"""

third_grid_voltage = clone_as_symbol(symbols.voltage, subscript="3")
"""
:symbols:`voltage` between the cathode and the third grid.
"""

anode_voltage = clone_as_symbol(symbols.voltage, display_symbol="U_a", display_latex="U_\\text{a}")
"""
:symbols:`voltage` between the cathode and the anode.
"""

first_grid_direct_permeability_coefficient = clone_as_symbol(symbols.direct_permeability_coefficient, subscript="1")
"""
:symbols:`direct_permeability_coefficient` of the first grid.
"""

second_grid_direct_permeability_coefficient = clone_as_symbol(symbols.direct_permeability_coefficient, subscript="2")
"""
:symbols:`direct_permeability_coefficient` of the second grid.
"""

third_grid_direct_permeability_coefficient = clone_as_symbol(symbols.direct_permeability_coefficient, subscript="3")
"""
:symbols:`direct_permeability_coefficient` of the third grid.
"""

anode_distance = clone_as_symbol(symbols.euclidean_distance, display_symbol="d_a", display_latex="d_\\text{a}")
"""
:symbols:`euclidean_distance` between the cathode and the anode.
"""

first_grid_distance = clone_as_symbol(symbols.euclidean_distance, subscript="1")
"""
:symbols:`euclidean_distance` between the cathode and the first grid.
"""

_first_expression = (
    first_grid_voltage
    + first_grid_direct_permeability_coefficient * second_grid_voltage
    + first_grid_direct_permeability_coefficient * second_grid_direct_permeability_coefficient * third_grid_voltage
)

_second_expression = (
    first_grid_direct_permeability_coefficient
    * second_grid_direct_permeability_coefficient
    * third_grid_direct_permeability_coefficient
    * anode_voltage
)

_third_expression = (
    1
    + ((anode_distance / first_grid_distance)**Rational(4, 3))
    * first_grid_direct_permeability_coefficient
)

law = Eq(equivalent_diode_voltage, (_first_expression + _second_expression) / _third_expression)
"""
:laws:symbol::

:laws:latex::
"""

@validate_input(voltage_of_first_grid_=first_grid_voltage,
    second_grid_voltage_=second_grid_voltage,
    third_grid_voltage_=third_grid_voltage,
    anode_voltage_=anode_voltage,
    coefficient_direct_permeability_of_first_grid_=first_grid_direct_permeability_coefficient,
    coefficient_direct_permeability_of_second_grid_=second_grid_direct_permeability_coefficient,
    coefficient_direct_permeability_of_third_grid_=third_grid_direct_permeability_coefficient,
    distance_to_anode_=anode_distance,
    distance_to_first_grid_=first_grid_distance)
@validate_output(equivalent_diode_voltage)
def calculate_voltage_of_equivalent_diode(voltage_of_first_grid_: Quantity,
    second_grid_voltage_: Quantity, third_grid_voltage_: Quantity, anode_voltage_: Quantity,
    coefficient_direct_permeability_of_first_grid_: float,
    coefficient_direct_permeability_of_second_grid_: float,
    coefficient_direct_permeability_of_third_grid_: float, distance_to_anode_: Quantity,
    distance_to_first_grid_: Quantity) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    if distance_to_anode_.scale_factor <= distance_to_first_grid_.scale_factor:
        raise ValueError(
            "The distance to the anode should be greater than the distance to the grid.")
    result_expr = solve(law, equivalent_diode_voltage, dict=True)[0][equivalent_diode_voltage]
    result_expr = result_expr.subs({
        first_grid_voltage:
            voltage_of_first_grid_,
        second_grid_voltage:
            second_grid_voltage_,
        third_grid_voltage:
            third_grid_voltage_,
        anode_voltage:
            anode_voltage_,
        first_grid_direct_permeability_coefficient:
            coefficient_direct_permeability_of_first_grid_,
        second_grid_direct_permeability_coefficient:
            coefficient_direct_permeability_of_second_grid_,
        third_grid_direct_permeability_coefficient:
            coefficient_direct_permeability_of_third_grid_,
        anode_distance:
            distance_to_anode_,
        first_grid_distance:
            distance_to_first_grid_
    })
    return Quantity(result_expr)
