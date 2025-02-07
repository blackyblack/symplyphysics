"""
Resonant frequency of ring resonator
====================================

The ring resonator is a microstrip line in the shape of a circle. When a wave propagates
along a microstrip line, part of the field goes out, since the microstrip line does not
have metal borders on all sides, unlike, for example, rectangular waveguides. A wave
traveling through an ring resonator acquires a phase shift and interacts with a wave
incident on the resonator. If the phase shift is expressed as :math:`2 \\pi N`, then
these waves add up in phase; here :math:`N` is the interference order.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities

frequency = symbols.temporal_frequency
"""
:symbols:`temporal_frequency` of the ring resonator.
"""

length = symbols.length
"""
:symbols:`length` of the circumference of the ring resonator.
"""

interference_order = symbols.positive_number
"""
Interference order, see :symbols:`positive_number`.
"""

relative_permittivity = symbols.relative_permittivity
"""
Effective :symbols:`relative_permittivity` of the resonator.
"""

law = Eq(frequency, interference_order * quantities.speed_of_light / (length * sqrt(relative_permittivity)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(ring_length_=length,
    order_interference_=interference_order,
    permittivity_=relative_permittivity)
@validate_output(frequency)
def calculate_frequency(ring_length_: Quantity, order_interference_: int,
    permittivity_: float) -> Quantity:
    result_expr = solve(law, frequency, dict=True)[0][frequency]
    result_expr = result_expr.subs({
        length: ring_length_,
        interference_order: order_interference_,
        relative_permittivity: permittivity_,
    })
    return Quantity(result_expr)
