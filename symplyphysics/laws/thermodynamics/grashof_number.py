from sympy import Eq, solve, S
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output, dimensionless, convert_to)


# Description
# The Grashof number (Gr) is a dimensionless number which the ratio of
# the buoyancy to viscous forces acting on a fluid. There is a characteristic
# length in the formula. The characteristic length is the dimension
# that defines the length scale of a physical system. A characteristic length
# is usually the volume of a system divided by its surface: L = V / A,
# where V is the volume of the body, and A is the cross-sectional area.
# For example, it is used to calculate flow through circular and non-circular
# tubes in order to examine flow conditions. D = 4 * A / p, where
# D is characteristic diameter, A is the cross-sectional are, p is wetted perimeter.
# Law: Gr = g * beta * (t_s - t_f) * L**3 / (nu ** 2), where
# g is gravitational acceleration,
# beta is coefficient of volume expansion,
# t_s is surface temperature,
# t_f is bulk temperature,
# L is characteristic length,
# nu is kinematic viscosity,
# Gr is Grashof number.

gravitational_acceleration = Symbol("acceleration_due_to_gravity", units.acceleration)
coefficient_of_volume_expansion = Symbol("coefficient_of_volume_expansion", units.temperature**-1)
surface_temperature = Symbol("initial_temperature", units.temperature)
bulk_temperature = Symbol("final_temperature", units.temperature)
characteristic_length = Symbol("characteristic_length", units.length)
viscosity = Symbol("viscosity", units.area / units.time)
grashof_number = Symbol("grashof_number", dimensionless)


law = Eq(grashof_number, gravitational_acceleration * coefficient_of_volume_expansion *
         (surface_temperature - bulk_temperature) * characteristic_length**3 / (viscosity**2))

def print_law() -> str:
    return print_expression(law)


@validate_input(
    acceleration_due_to_gravity_=gravitational_acceleration,
    coefficient_of_volume_expansion_=coefficient_of_volume_expansion,
    surface_temperature_=surface_temperature,
    bulk_temperature_=bulk_temperature,
    characteristic_length_=characteristic_length,
    viscosity_=viscosity
)
@validate_output(grashof_number)
def calculate_grashof_number(
    acceleration_due_to_gravity_: Quantity,
    coefficient_of_volume_expansion_: Quantity,
    surface_temperature_: Quantity,
    bulk_temperature_: Quantity,
    characteristic_length_: Quantity,
    viscosity_: Quantity
) -> float:
    result_expr = solve(law, grashof_number, dict=True)[0][grashof_number]
    result_applied = result_expr.subs({
        gravitational_acceleration: acceleration_due_to_gravity_,
        coefficient_of_volume_expansion: coefficient_of_volume_expansion_,
        surface_temperature: surface_temperature_,
        bulk_temperature: bulk_temperature_,
        characteristic_length: characteristic_length_,
        viscosity: viscosity_
    })
    result = Quantity(result_applied)
    return float(convert_to(result, S.One).evalf())