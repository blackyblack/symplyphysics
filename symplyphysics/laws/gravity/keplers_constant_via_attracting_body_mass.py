"""
Kepler's constant via attracting body mass
==========================================

*Kepler's constant* is a physical quantity that only depends on the gravitational constant and the mass of
the orbited body, such as the Sun. It is constant in the sense that all planets that orbit the same body
approximately have the same value of the Kepler's constant.
"""

from sympy import Eq, pi
from sympy.physics.units import gravitational_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

keplers_constant = Symbol("keplers_constant", units.length**3 / units.time**2)
r"""
Symbol:
    K

Latex:
    :math:`\mathfrak{K}`
"""

attracting_mass = clone_symbol(symbols.basic.mass, display_symbol="M")
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the attracting body.
"""

law = Eq(keplers_constant, gravitational_constant * attracting_mass / (4 * pi**2))
r"""
K = G * M / (4 * pi^2)

Latex:
    :math:`\mathfrak{K} = \frac{G M}{4 \pi^2}`
"""


@validate_input(attracting_body_mass_=attracting_mass)
@validate_output(keplers_constant)
def calculate_keplers_constant(attracting_body_mass_: Quantity) -> Quantity:
    result = law.rhs.subs({
        attracting_mass: attracting_body_mass_,
    })
    return Quantity(result)
