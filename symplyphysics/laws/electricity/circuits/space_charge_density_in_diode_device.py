from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)
from sympy.physics.units import electric_constant

# Description
## The simplest device implementing a cathode current control method is a diode device.
## It contains a thermionic cathode, a mesh electrode and an anode.

## Law is: p = (4 / 9) * eps0 * U / d^2, where
## p - space charge density in the diode device,
## eps0 - dielectric constant,
## U - voltage on the grid,
## d - distance between the grid and the cathode.

space_charge_density = Symbol("space_charge_density", units.charge / units.volume)

voltage_on_grid = Symbol("voltage_on_grid", units.voltage)
distance_between_cathode_and_grid = Symbol("distance_between_cathode_and_grid", units.length)

law = Eq(space_charge_density,
    (4 / 9) * electric_constant * voltage_on_grid / distance_between_cathode_and_grid**2)


@validate_input(voltage_on_grid_=voltage_on_grid,
    distance_between_cathode_and_grid_=distance_between_cathode_and_grid)
@validate_output(space_charge_density)
def calculate_space_charge_density(voltage_on_grid_: Quantity,
    distance_between_cathode_and_grid_: Quantity) -> Quantity:
    result_expr = solve(law, space_charge_density, dict=True)[0][space_charge_density]
    result_expr = result_expr.subs({
        voltage_on_grid: voltage_on_grid_,
        distance_between_cathode_and_grid: distance_between_cathode_and_grid_,
    })
    return Quantity(result_expr)
