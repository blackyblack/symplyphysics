"""
Power carried by coaxial waveguide
==================================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt, ln
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

power = symbols.power
"""
:symbols:`power` transmitted by the waveguide.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the insulator.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the insulator.
"""

voltage = symbols.voltage
"""
:symbols:`voltage` between the central conductor and the outer conductor.
"""

outer_diameter = clone_as_symbol(symbols.diameter, display_symbol="d_o", display_latex="d_\\text{o}")
"""
:symbols:`diameter` of the outer conductor.
"""

inner_diameter = clone_as_symbol(symbols.diameter, display_symbol="d_i", display_latex="d_\\text{i}")
"""
:symbols:`diameter` of the inner conductor.
"""

impedance_constant = Quantity(120 * units.ohm, display_symbol="Z_0")
"""
Constant equal to :math:`120 \\Omega`.

..
    rename back to `vacuum_impedance`?
"""

law = Eq(power, (voltage**2 / impedance_constant) * sqrt(relative_permittivity /
    (relative_permeability * ln(outer_diameter / inner_diameter))))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability,
    voltage_=voltage,
    outer_diameter_=outer_diameter,
    inner_diameter_=inner_diameter)
@validate_output(power)
def calculate_waveguide_power(relative_permittivity_: float, relative_permeability_: float,
    voltage_: Quantity, outer_diameter_: Quantity, inner_diameter_: Quantity) -> Quantity:
    if outer_diameter_.scale_factor <= inner_diameter_.scale_factor:
        raise ValueError("The outer diameter must be greater than the inner diameter")
    result_expr = solve(law, power, dict=True)[0][power]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
        voltage: voltage_,
        outer_diameter: outer_diameter_,
        inner_diameter: inner_diameter_
    })
    return Quantity(result_expr)
