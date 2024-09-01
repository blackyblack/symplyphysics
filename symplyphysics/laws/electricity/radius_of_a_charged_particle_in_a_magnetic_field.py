from sympy import (Eq, solve, pi)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import period_of_a_charged_particle_in_a_magnetic_field as period_law
from symplyphysics.laws.kinematics import position_via_constant_speed_and_time as distance_law

# Description
## Let an arbitrary particle move in a magnetic field around a circle. Then radius of its motion depends on
## magnetic induction and on velocity, charge and mass of the particle.

## Law is: r = m * v / (q * B), where
## r - radius,
## m - mass,
## v - velocity,
## q - charge,
## B - induction.

radius = Symbol("radius", units.length)

velocity = Symbol("velocity", units.velocity)
charge = Symbol("charge", units.charge)
induction = Symbol("induction", units.magnetic_density)
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass")

law = Eq(radius, particle_mass * velocity / (charge * induction))

# This law might be derived via period of a charged particle in a magnetic field and distance from constant velocity.

distance_law_applied = distance_law.law.subs({
    distance_law.initial_position: 0,
    distance_law.final_position: 2 * pi * radius,
    distance_law.speed: velocity,
})
# Period is a time taken to cover a full circle.
period_derived = solve(distance_law_applied, distance_law.time,
    dict=True)[0][distance_law.time]

law_applied = period_law.law.subs({
    period_law.particle_mass: particle_mass,
    period_law.charge: charge,
    period_law.induction: induction,
    period_law.period: period_derived,
})
radius_derived = solve(law_applied, radius, dict=True)[0][radius]

# Check if derived radius is same as declared.
assert expr_equals(radius_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=particle_mass, velocity_=velocity, induction_=induction, charge_=charge)
@validate_output(radius)
def calculate_radius(mass_: Quantity, velocity_: Quantity, induction_: Quantity,
    charge_: Quantity) -> Quantity:
    result_expr = solve(law, radius, dict=True)[0][radius]
    result_expr = result_expr.subs({
        particle_mass: mass_,
        velocity: velocity_,
        induction: induction_,
        charge: charge_,
    })
    return Quantity(result_expr)
