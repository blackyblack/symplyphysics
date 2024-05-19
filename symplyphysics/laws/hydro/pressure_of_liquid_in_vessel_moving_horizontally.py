from sympy import (Eq, solve, sqrt)
from sympy.physics.units import acceleration_due_to_gravity as earth_free_fall_acceleration
from symplyphysics import (
    Vector,
    add_cartesian_vectors,
    symbols,
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    vector_magnitude,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.hydro import hydrostatic_pressure_from_density_and_depth_acceleration as pressure_law

# Description
## If a vessel with a liquid moves with horizontal acceleration, then the hydrostatic pressure of the liquid
## depends on the density of the liquid, the acceleration of free fall, the horizontal acceleration of the vessel
## and the height of liquid.
## Height is the distance from the surface of the water in the vessel to the surface of equal pressure. In a vessel
## moving with horizontal acceleration, the planes of equal pressure are inclined at an angle to the horizon.

## Law is: p = p0 * sqrt(g^2 + a^2) * h, where
## p - pressure,
## p0 - density of liquid,
## g - earth free fall acceleration,
## a - acceleration of vessel,
## h - height.

pressure = Symbol("pressure", units.pressure)

density_liquid = Symbol("density_liquid", units.mass / units.volume)
acceleration = symbols.kinematic.acceleration
height = Symbol("height", units.length)

law = Eq(pressure,
    density_liquid * sqrt(earth_free_fall_acceleration**2 + acceleration**2) * height)

# This law might be derived via hydrostatic pressure law.
# The vessel moves horizontally and the pressure exerted by the resultant force on a surface of equal pressure,
# which is inclined with respect to the horizon, is considered. The modulus of the optimizing force will be equal
# to m * sqrt(g^2+ a^2).

# Earth free fall acceleration is directed downwards. It is equivalent to vessel acceleration
# directed vertically upwards.
free_fall_acceleration_vector = Vector([0, earth_free_fall_acceleration])
# Horizontal vector
vessel_acceleration_vector = Vector([acceleration, 0])
total_acceleration = vector_magnitude(
    add_cartesian_vectors(free_fall_acceleration_vector, vessel_acceleration_vector))

pressure_law_applied = pressure_law.law.subs({
    pressure_law.density: density_liquid,
    pressure_law.depth: height,
    pressure_law.acceleration: total_acceleration,
})
pressure_derived = solve(pressure_law_applied, pressure_law.hydrostatic_pressure,
    dict=True)[0][pressure_law.hydrostatic_pressure]

# Check if derived pressure is same as declared.
assert expr_equals(pressure_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(density_liquid_=density_liquid, acceleration_=acceleration, height_=height)
@validate_output(pressure)
def calculate_pressure(density_liquid_: Quantity, acceleration_: Quantity,
    height_: Quantity) -> Quantity:
    result_expr = solve(law, pressure, dict=True)[0][pressure]
    result_expr = result_expr.subs({
        density_liquid: density_liquid_,
        acceleration: acceleration_,
        height: height_
    })
    return Quantity(result_expr)
