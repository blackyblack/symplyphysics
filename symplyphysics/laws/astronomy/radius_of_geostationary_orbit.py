from sympy import (Eq, Rational, solve)
from sympy.physics.units import gravitational_constant
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)

# Description
## A geostationary orbit is a circular orbit located above the Earth's equator (0° latitude), where an artificial satellite orbits the planet
## with an angular velocity equal to the angular velocity of the Earth's rotation around its axis.

## Law is: R = (G * M / (w^2))^(1/3), where
## R - radius of geostationary orbit,
## G - gravitational constant,
## M - mass of planet.
## w - angular speed of rotation of satellite.

# Links: possible formula derivable from here <https://en.wikipedia.org/wiki/Geostationary_orbit#Derivation>
# TODO: find link with exact formula

# TODO: move law to `gravity`?

radius_of_orbit = Symbol("radius_of_orbit", units.length)
mass_of_planet = symbols.mass
speed_rotation_satellite = Symbol("speed_rotation_satellite", angle_type / units.time)

law = Eq(radius_of_orbit,
    (gravitational_constant * mass_of_planet / (speed_rotation_satellite**2))**Rational(1, 3))


@validate_input(mass_of_planet_=mass_of_planet, speed_rotation_satellite_=speed_rotation_satellite)
@validate_output(radius_of_orbit)
def calculate_radius_of_orbit(mass_of_planet_: Quantity,
    speed_rotation_satellite_: Quantity) -> Quantity:
    result_expr = solve(law, radius_of_orbit, dict=True)[0][radius_of_orbit]
    result_expr = result_expr.subs({
        mass_of_planet: mass_of_planet_,
        speed_rotation_satellite: speed_rotation_satellite_,
    })
    return Quantity(result_expr)
