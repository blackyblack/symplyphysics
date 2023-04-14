from symplyphysics import (
    expr_to_quantity, units, convert_to,solve, SI, Eq, pretty
)

from symplyphysics.laws.electricity import amount_energy_from_voltage_time_resistance as joule_lenz_law
from symplyphysics.laws.thermodynamics import thermal_energy_from_mass_and_temperature as operate_energy
from symplyphysics.definitions import density_from_mass_volume as density_law
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohm_law
from symplyphysics.definitions import power_is_proportional_voltage_and_current as operate_power

# The household electric kettle heated 0.5 liters of water from 20 degree Celsius to boiling.
# The power of the kettle is 1500 watts. How long did the heating process take?
# Since physics uses Kelvin temperature units, we must consider that 0 degrees Celsius = 273 Kelvin.

CELSIUS_TO_KELVIN_OFFSET = int(273)

# water parameters: density and specific heat capacity
water_density = units.Quantity('water_density')
SI.set_quantity_dimension(water_density, units.mass / units.volume)
SI.set_quantity_scale_factor(water_density, 1000 * units.kilogram / units.meter**3)
water_heat_capacity = units.Quantity('water_heat_capacity')
SI.set_quantity_dimension(water_heat_capacity, units.energy / (units.mass * units.temperature))
SI.set_quantity_scale_factor(water_heat_capacity, 4200 * units.joule / (units.kilogram * units.kelvin))

# kettle parameters: power and volume
kettle_power = units.Quantity('kettle_power')
SI.set_quantity_dimension(kettle_power, units.power)
SI.set_quantity_scale_factor(kettle_power, 1500 * units.watt)
kettle_volume = units.Quantity('kettle_volume')
SI.set_quantity_dimension(kettle_volume, units.volume)
SI.set_quantity_scale_factor(kettle_volume, 0.5 * units.liter)

# heating parameters
initial_temperature = units.Quantity('initial_temperature')
SI.set_quantity_dimension(initial_temperature, units.temperature)
SI.set_quantity_scale_factor(initial_temperature, (20 + CELSIUS_TO_KELVIN_OFFSET) * units.kelvin)
final_temperature = units.Quantity('final_temperature')
SI.set_quantity_dimension(final_temperature, units.temperature)
SI.set_quantity_scale_factor(final_temperature, (100 + CELSIUS_TO_KELVIN_OFFSET) * units.kelvin)

ohm_law_applied = ohm_law.law.subs({ohm_law.current: operate_power.current, ohm_law.resistance: joule_lenz_law.resistance, ohm_law.voltage: joule_lenz_law.voltage})
operate_power_applied = operate_power.law.subs(operate_power.voltage, joule_lenz_law.voltage)
density_applied = density_law.definition.subs(density_law.mass, operate_energy.body_mass)

law = [density_applied, ohm_law_applied, operate_power_applied, operate_energy.law, joule_lenz_law.law]
heating_time_solved = solve(law, (joule_lenz_law.amount_energy, joule_lenz_law.resistance, operate_power.current, operate_energy.body_mass, joule_lenz_law.time), dict=True)[0][joule_lenz_law.time]
heating_time_eq = Eq(joule_lenz_law.time, heating_time_solved)

print("\nFormula is:\n\n {}".format(pretty(heating_time_eq, use_unicode=False)))

heating_time_expr = heating_time_solved.subs({
    operate_power.power: kettle_power,
    density_law.density: water_density,
    density_law.volume: kettle_volume,
    operate_energy.specific_heat_capacity: water_heat_capacity,
    operate_energy.temperature_origin: initial_temperature, operate_energy.temperature_end: final_temperature
})
heating_time = expr_to_quantity(heating_time_expr, "heating_time")

kettle_power_value = convert_to(kettle_power, units.watt).subs(units.watt, 1).evalf(5)
kettle_volume_value = convert_to(kettle_volume, units.liter).subs(units.liter, 1).evalf(3)
kettle_t1_value = convert_to(initial_temperature, units.kelvin).subs(units.kelvin, 1).evalf(3) - CELSIUS_TO_KELVIN_OFFSET
kettle_t2_value = convert_to(final_temperature, units.kelvin).subs(units.kelvin, 1).evalf(3) - CELSIUS_TO_KELVIN_OFFSET
heating_time_value = convert_to(heating_time, units.second).subs(units.second, 1).evalf(3)

print(f"\nPower = {kettle_power_value} {units.watt}, water volume = {kettle_volume_value} {units.liter}, temperature_begin = {kettle_t1_value} celsius degrees, temperature_end = {kettle_t2_value} celsius degrees: heating time = {heating_time_value} {units.second}\n")
print(f"\nSolution: heating time = {heating_time_value} {units.second}\n")
