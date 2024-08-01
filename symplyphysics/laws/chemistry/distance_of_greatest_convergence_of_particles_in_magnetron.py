from sympy import Eq, Rational, solve, log
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
)

# Description
## The traveling atom moves towards the substrate in the magnetron. At the same time, it collides with gas atoms.
## The energy transfer coefficient in these collisions depends on the mass of the traveling atom and the mass of the gas atom.
## The distance of the greatest convergence of colliding particles can be calculated using the model of quasi-rigid spheres.
## The discharge voltage is the voltage between the cathode and the anode in the magnetron at which plasma occurs.

## Law is: r = -(0.122e-10 * (Z1^0.0387 + Z2^0.0387)) * ln(U / (95.863 * (Z1 * Z2)^0.7383 * (Z1 * Z2)^0.75), where
## r - the distance of the greatest convergence of two colliding particles,
## U - discharge voltage,
## Z1 - the atomic number of the first atom (traveling atom),
## Z2 - the atomic number of the second atom (gas atom).

distance_of_convergence_of_particles = Symbol("distance_of_convergence_of_particles", units.length)

discharge_voltage = Symbol("discharge_voltage", units.voltage)
atomic_number_of_first_atom = Symbol("atomic_number_of_first_atom", dimensionless)
atomic_number_of_second_atom = Symbol("atomic_number_of_first_atom", dimensionless)

first_constant = Quantity(0.122e-10 * units.meter)
second_constant = Quantity(95.863 * units.volt)
expression_1 = first_constant * (atomic_number_of_first_atom**Rational(0.0387) +
    atomic_number_of_second_atom**Rational(0.0387))
expression_2 = second_constant * (atomic_number_of_first_atom *
    atomic_number_of_second_atom)**Rational(0.7383)

law = Eq(
    distance_of_convergence_of_particles, -expression_1 * log(discharge_voltage / (expression_2 *
    (atomic_number_of_first_atom * atomic_number_of_second_atom)**Rational(0.75))))


@validate_input(discharge_voltage_=discharge_voltage,
    atomic_number_of_first_atom_=atomic_number_of_first_atom,
    atomic_number_of_second_atom_=atomic_number_of_second_atom)
@validate_output(distance_of_convergence_of_particles)
def calculate_distance_of_convergence_of_particles(discharge_voltage_: Quantity,
    atomic_number_of_first_atom_: int, atomic_number_of_second_atom_: int) -> Quantity:
    result_expr = solve(law, distance_of_convergence_of_particles,
        dict=True)[0][distance_of_convergence_of_particles]
    result_expr = result_expr.subs({
        discharge_voltage: discharge_voltage_,
        atomic_number_of_first_atom: atomic_number_of_first_atom_,
        atomic_number_of_second_atom: atomic_number_of_second_atom_,
    })
    return Quantity(result_expr)
