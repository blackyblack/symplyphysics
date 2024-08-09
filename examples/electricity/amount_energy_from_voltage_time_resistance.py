#!/usr/bin/env python3

from sympy import solve, Eq, dsolve
from symplyphysics import (units, convert_to, Quantity, prefixes)
from symplyphysics.core.symbols.celsius import Celsius, to_kelvin_quantity
from symplyphysics.definitions import (
    density_from_mass_volume as density_def,
    power_is_energy_derivative as power_def,
)
from symplyphysics.laws.thermodynamics import (
    heat_is_heat_capacity_times_temperature_change as thermal_energy_law,)
from symplyphysics.laws.quantities import quantity_is_specific_quantity_times_mass as specific_qty_law

# The household electric kettle heated 0.5 liters of water from 20 degree Celsius to boiling.
# The power of the kettle is 1500 watts. How long did the heating process take?
# Since physics uses Kelvin temperature units, we must convert Celsius degrees to Kelvin.

# water parameters: density and specific heat capacity
water_density = Quantity(1000 * units.kilogram / units.meter**3)
water_heat_capacity = Quantity(4.2 * prefixes.kilo * units.joule / (units.kilogram * units.kelvin))

# kettle parameters: power and volume
kettle_power = Quantity(1500 * units.watt)
kettle_volume = Quantity(0.5 * units.liter)

# heating parameters
initial_temperature = Celsius(20)
final_temperature = Celsius(100)
initial_temperature_kelvin = to_kelvin_quantity(initial_temperature)
final_temperature_kelvin = to_kelvin_quantity(final_temperature)

water_mass = solve(density_def.definition, density_def.mass)[0].subs({
    density_def.density: water_density,
    density_def.volume: kettle_volume,
})

heat_via_power = dsolve(
    power_def.definition.subs(power_def.power(power_def.time), kettle_power),
    power_def.energy(power_def.time),
    ics={power_def.energy(0): 0},
).rhs

water_heat_capacity = specific_qty_law.law.rhs.subs({
    specific_qty_law.specific_quantity: water_heat_capacity,
    specific_qty_law.mass: water_mass,
})

water_temperature_change = final_temperature_kelvin - initial_temperature_kelvin

heat_transferred = thermal_energy_law.law.rhs.subs({
    thermal_energy_law.heat_capacity: water_heat_capacity,
    thermal_energy_law.temperature_change: water_temperature_change,
})

heating_time = solve(
    Eq(heat_via_power, heat_transferred),
    power_def.time,
)[0]

heating_time_value = convert_to(Quantity(heating_time), units.second).evalf(3)

print(f"Heating time is {heating_time_value} s.")
