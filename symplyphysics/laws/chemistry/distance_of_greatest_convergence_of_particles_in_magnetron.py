from sympy import Eq, solve, log
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, dimensionless,)

# Description
## The traveling atom moves towards the substrate in the magnetron. At the same time, it collides with gas atoms.
## The energy transfer coefficient in these collisions depends on the mass of the traveling atom and the mass of the gas atom.
## The distance of the greatest convergence of colliding particles can be calculated using the model of quasi-rigid spheres.

## Law is: r = -(0.122e-10 * 2 * Z^0.0387) * ln(U / (95.863 * Z^1.4766 * (Z^2)^0.75), where
## r - the distance of the greatest convergence of two colliding particles,
## U - discharge voltage,
## Z - the number of the gas atom.

distance_of_convergence_of_particles = Symbol("distance_of_convergence_of_particles", units.length)

discharge_voltage = Symbol("discharge_voltage", units.voltage)
number_of_gas_atom = Symbol("number_of_gas_atom", dimensionless)

first_constant = Quantity(0.122e-10 * units.meter)
second_constant = Quantity(95.863 * units.volt)
expression_1 = first_constant * 2 * number_of_gas_atom**0.0387
expression_2 = second_constant * number_of_gas_atom**1.4766

law = Eq(distance_of_convergence_of_particles, -expression_1 * log(discharge_voltage / (expression_2 * (number_of_gas_atom**2)**0.75)))


@validate_input(discharge_voltage_=discharge_voltage,
    number_of_gas_atom_=number_of_gas_atom)
@validate_output(distance_of_convergence_of_particles)
def calculate_distance_of_convergence_of_particles(discharge_voltage_: Quantity,
    number_of_gas_atom_: float) -> Quantity:
    result_expr = solve(law, distance_of_convergence_of_particles, dict=True)[0][distance_of_convergence_of_particles]
    result_expr = result_expr.subs({
        discharge_voltage: discharge_voltage_,
        number_of_gas_atom: number_of_gas_atom_,
    })
    return Quantity(result_expr)
