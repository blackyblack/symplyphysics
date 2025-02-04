"""
Current from voltage and triode constant in vacuum triode
=========================================================

A triode has three electrodes: a cathode, an anode and one control grid. The triode can
be replaced with an equivalent diode and the :math:`3/2`-power law can be applied. The
triode constant in this law depends only on the relative position, shape and size of the
electrodes of the vacuum triode.

**Links**

#. `StudFiles, formula after (1.15) <https://studfile.net/preview/7192538/page:5/>`__.

..
    TODO: find link in English
"""

from sympy import Eq, Rational, solve
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

anode_current = symbols.current
"""
Anode :symbols:`current`.
"""

triode_constant = SymbolNew("g", units.current / units.voltage**Rational(3, 2))
"""
Triode constant which depends only on the relative position, shape and size of the
electrodes of the vacuum triode.
"""

anode_voltage = clone_as_symbol(symbols.voltage, display_symbol="U_a", display_latex="U_\\text{a}")
"""
:symbols:`voltage` between cathode and anode.
"""

voltage_gain = clone_as_symbol(symbols.circuit_gain, subscript="V")
"""
Voltage gain in the triode. See :symbols:`circuit_gain`.
"""

grid_voltage = clone_as_symbol(symbols.voltage, display_symbol="U_g", display_latex="U_\\text{g}")
"""
:symbols:`voltage` between cathode and grid.
"""

law = Eq(anode_current,
    triode_constant * (anode_voltage + voltage_gain * grid_voltage)**Rational(3, 2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(triode_constant_=triode_constant,
    anode_voltage_=anode_voltage,
    voltage_triode_gain_=voltage_gain,
    grid_voltage_=grid_voltage)
@validate_output(anode_current)
def calculate_current(triode_constant_: Quantity, anode_voltage_: Quantity,
    voltage_triode_gain_: float, grid_voltage_: Quantity) -> Quantity:
    result_expr = solve(law, anode_current, dict=True)[0][anode_current]
    result_expr = result_expr.subs({
        triode_constant: triode_constant_,
        anode_voltage: anode_voltage_,
        voltage_gain: voltage_triode_gain_,
        grid_voltage: grid_voltage_,
    })
    return Quantity(result_expr)
