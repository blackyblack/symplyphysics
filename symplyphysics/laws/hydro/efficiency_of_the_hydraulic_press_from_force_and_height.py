"""
Efficiency of hydraulic press from force and height
===================================================

A real hydraulic press is never :math:`100%` efficient due to friction and other energy
losses. Its efficiency is the ratio of the useful work (given by the product of output
force and the height of the output side) to the expended work (given by the product of
input force and the height of the input side).

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
)

output_force = clone_as_symbol(symbols.force, subscript="2")
"""
Output :symbols:`force`.
"""

output_distance = clone_as_symbol(symbols.euclidean_distance, subscript="2")
"""
:symbols:`euclidean_distance` covered by the output piston.
"""

input_force = clone_as_symbol(symbols.force, subscript="1")
"""
Input :symbols:`force`.
"""

input_distance = clone_as_symbol(symbols.euclidean_distance, subscript="1")
"""
:symbols:`euclidean_distance` covered by the input piston.
"""

efficiency = symbols.mechanical_efficiency
"""
:symbols:`mechanical_efficiency` of the hydraulic press.
"""

law = Eq(efficiency, (output_force * output_distance) / (input_force * input_distance))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(useful_force_=output_force,
    useful_height_=output_distance,
    expended_force_=input_force,
    expended_height_=input_distance)
@validate_output(efficiency)
def calculate_efficiency(useful_force_: Quantity, useful_height_: Quantity,
    expended_force_: Quantity, expended_height_: Quantity) -> float:
    result_expr = solve(law, efficiency, dict=True)[0][efficiency]
    result_efficiency = result_expr.subs({
        output_force: useful_force_,
        output_distance: useful_height_,
        input_force: expended_force_,
        input_distance: expended_height_
    })
    return convert_to_float(result_efficiency)
