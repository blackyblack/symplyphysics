from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

# Description
## The law of areas, also known as Kepler's second law of planetary motion, states that
## a line the connects a planet to the attracting body (the Sun) sweeps out equal areas
## in the plane of the planet's orbit in equal time intervals, and its rate of change
## is proportional to the planets' angular momentum.

# Note: this law is equivalent to saying that the planet's angular momentum is conserved.

# Law: dA/dt = L/(2*m)
## A - sweeping area
## d/dt - derivative w.r.t. time
## L - planet's angular momentum
## m - planet's mass

time = Symbol("time", units.time)
area_swept = Function("area_swept", units.area)
planet_angular_momentum = Symbol("planet_angular_momentum", units.length * units.momentum)
planet_mass = symbols.mass

law = Eq(Derivative(area_swept(time), time), planet_angular_momentum / (2 * planet_mass))

# TODO: derive law from expression of area and definition of angular momentum


def print_law() -> str:
    return print_expression(law)


@validate_input(planet_angular_momentum_=planet_angular_momentum, planet_mass_=planet_mass)
@validate_output(units.area / units.time)
def calculate_rate_of_change_of_area(planet_angular_momentum_: Quantity,
    planet_mass_: Quantity) -> Quantity:
    result = law.rhs.subs({
        planet_angular_momentum: planet_angular_momentum_,
        planet_mass: planet_mass_,
    })
    return Quantity(result)
