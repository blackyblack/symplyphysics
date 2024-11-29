"""
Area rate of change is proportional to angular momentum
=======================================================

The **law of areas**, also known as **Kepler's second law of planetary motion**,
states that a line the connects a planet to the attracting body (the Sun) sweeps
out equal areas in the plane of the planet's orbit in equal time intervals, and
its rate of change is proportional to the planet's angular momentum. It is
equivalent to saying that the planet's angular momentum is conserved.

#. Sivukhin D.V. (1979), __Obshchiy kurs fiziki__ [General course of Physics], vol. 1, pp. 312—314.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)

time = symbols.time
"""
:symbols:`time`.
"""

area_swept = clone_as_function(symbols.area, [time])
"""
:symbols:`area` swept by the planet.
"""

planet_angular_momentum = symbols.angular_momentum
"""
:symbols:`angular_momentum` of the planet.
"""

planet_mass = symbols.mass
"""
:symbols:`mass` of the planet.
"""

law = Eq(Derivative(area_swept(time), time), planet_angular_momentum / (2 * planet_mass))
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive law from expression of area and definition of angular momentum


@validate_input(planet_angular_momentum_=planet_angular_momentum, planet_mass_=planet_mass)
@validate_output(units.area / units.time)
def calculate_rate_of_change_of_area(planet_angular_momentum_: Quantity,
    planet_mass_: Quantity) -> Quantity:
    result = law.rhs.subs({
        planet_angular_momentum: planet_angular_momentum_,
        planet_mass: planet_mass_,
    })
    return Quantity(result)
