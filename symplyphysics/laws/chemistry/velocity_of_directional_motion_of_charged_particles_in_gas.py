from sympy import Eq, solve
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)

# Description
## In a gas, the free path length of charged particles is small and the particles, as they move, experience many collisions in the
## volume between the electrodes. In this case, the speed of directional movement will depend not on the potential difference traveled,
## but on the local value of the electric intensity.

## Law is: v = b0 * E / p, where
## v - the velocity of the directional motion of charged particles in a gas,
## b0 - mobility of charged particles at a single pressure,
## E - electric intensity,
## p - pressure.

velocity = Symbol("velocity", units.velocity)

mobility_at_unit_pressure = Symbol("mobility_at_unit_pressure", units.velocity * units.pressure * units.length / units.voltage)
pressure = Symbol("pressure", units.pressure)
electric_intensity = Symbol("pressure", units.voltage / units.length)

law = Eq(velocity, mobility_at_unit_pressure * electric_intensity / pressure)


@validate_input(mobility_at_unit_pressure_=mobility_at_unit_pressure,
    pressure_=pressure,
    electric_intensity_=electric_intensity)
@validate_output(velocity)
def calculate_velocity(mobility_at_unit_pressure_: Quantity,
    pressure_: Quantity, electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, velocity,
        dict=True)[0][velocity]
    result_expr = result_expr.subs({
        mobility_at_unit_pressure: mobility_at_unit_pressure_,
        pressure: pressure_,
        electric_intensity: electric_intensity_,
    })
    return Quantity(result_expr)
